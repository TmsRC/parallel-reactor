#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include "simulation_configuration.h"
#include "simulation_support.h"

#include <limits.h>
#include <mpi.h>

// Height of a fuel pellet in meters (they are 40mm by 40mm by 2.5mm)
#define HEIGHT_FUEL_PELLET_M 0.0025
// Size of each reactor channel in meters (they are cuboid, so this value in x and y)
#define CHANNEL_SIZE 0.2
// Depth of fuel pellets in meters (they are x=40mm by y=40mm by z=2mm)
#define FUEL_PELLET_DEPTH 0.002
// Weight in grams of the neutron generator for each cm in height
#define NEUTRON_GENERATOR_WEIGHT_PER_CM 0.5
// Number of nanoseconds in a second
#define NS_AS_SEC 1e-9

struct fission_event {
    unsigned long int generated_neutrons;
    int pellet_height;
};


// The neutrons that are currently moving throughout the reactor core
struct neutron_struct *neutrons;
// Indexes of empty neutrons (i.e. those inactive) that can be used
unsigned long int *neutron_index;
// The current index in the neutron_index array which is the next free neutron
unsigned long int currentNeutronIndex = 0;
// The reactor core itself, each are channels in the x and y dimensions
struct channel_struct **reactor_core;

/*static*/ void step(int, struct simulation_configuration_struct *);
/*static*/ void generateReport(int, int, struct simulation_configuration_struct *, struct timeval);
/*static*/ void updateReactorCore(int, struct simulation_configuration_struct *);
/*static*/ void updateNeutrons(int, struct simulation_configuration_struct *);
/*static*/ void updateFuelAssembly(int, struct channel_struct *);
/*static*/ void updateNeutronGenerator(int, struct channel_struct *, struct simulation_configuration_struct *);
/*static*/ void createNeutrons(int, struct channel_struct *, double);
/*static*/ void initialiseReactorCore(struct simulation_configuration_struct *);
/*static*/ void initialiseNeutrons(struct simulation_configuration_struct *);
/*static*/ double getControlRodLoweredToLevel(struct simulation_configuration_struct *, int, int);
/*static*/ void writeReactorState(struct simulation_configuration_struct *, int, char *);
/*static*/ void getFuelAssemblyChemicalContents(struct fuel_assembly_struct *, double *);
/*static*/ void clearReactorStateFile(char *);
/*static*/ struct channel_struct *locateChannelFromPosition(double, double, struct simulation_configuration_struct *);
/*static*/ unsigned long int getTotalNumberFissions(struct simulation_configuration_struct *);
/*static*/ unsigned long int getNumberActiveNeutrons(struct simulation_configuration_struct *);
/*static*/ double getElapsedTime(struct timeval);

void updateNeutronPosition(struct neutron_struct *neutron, int dt);
void interactWithFuelAssembly(struct neutron_struct *neutron, struct channel_struct *reactorChannel, struct simulation_configuration_struct *configuration, long int i);
void interactWithModerator(struct neutron_struct *neutron, struct channel_struct *reactorChannel, struct simulation_configuration_struct *configuration, long int i);
void interactWithControlRod(struct neutron_struct *neutron, struct channel_struct *reactorChannel, struct simulation_configuration_struct *configuration, long int i);
void interactWithReactor(struct neutron_struct *neutron, struct channel_struct *reactorChannel, struct simulation_configuration_struct *configuration, long int i);
bool determineOutbound(struct neutron_struct *neutron, struct simulation_configuration_struct *configuration, long int i);
void rebuildIndex(struct simulation_configuration_struct *configuration);

void commit_simple_neutron_datatype(void);
void commit_fission_event_datatype(void);
void commit_sparse_neutrons_datatype(int count);


void executeFissions(int dt, struct channel_struct *channel);
void manageFuelAssemblyFissions(int dt, struct simulation_configuration_struct *configuration);
void createNeutronsFromFission(struct channel_struct *channel);
unsigned long int calculateNumberActiveNeutrons(struct simulation_configuration_struct *);

// MPI variables
int size, rank;
MPI_Comm world;
MPI_Datatype MPI_SIMPLE_NEUTRON;
MPI_Datatype MPI_SPARSE_NEUTRONS;
MPI_Datatype MPI_FISSION_EVENT;

MPI_Request neutrons_send_request;
MPI_Request neutrons_recv_request;

// Temporary variables

/*long*/ int *fuel_assembly_neutrons_index; // Note: I have to address this datatype issue
int sender, receiver;
unsigned long int fission_neutrons_counter;
int max_pellets;
int max_fission_events;
struct fission_event *fission_array;







void manageFuelAssemblyInteractions(struct simulation_configuration_struct *configuration)
{
    MPI_Status status;
    int count;
    MPI_Irecv(&neutrons[0],configuration->max_neutrons,MPI_SIMPLE_NEUTRON,sender,0,world,&neutrons_recv_request); //Note: Need to change this so that I can use a simpler struct than neutrons in rank 0
    MPI_Wait(&neutrons_recv_request,&status);
//    printf("Rank %d -- Received neutrons.\n",rank);

    MPI_Get_count(&status,MPI_SIMPLE_NEUTRON,&count);

//    printf("Rank %d, count %ld\n",rank,count);

//    for(int i=0; i<count; i+=100) printf("Rank %d, neutron %d, %lf,\n",rank,i,neutrons[i].energy);

    for(long int i=0; i<count; i++)
    {
            struct channel_struct *reactorChannel = locateChannelFromPosition(neutrons[i].pos_x, neutrons[i].pos_y, configuration);

            interactWithFuelAssembly(&neutrons[i],reactorChannel,configuration,0);

    }

//    printf("Rank %d -- Sending back...\n",rank);
    MPI_Issend(&neutrons[0],count,MPI_SIMPLE_NEUTRON,sender,0,world,&neutrons_send_request);
    MPI_Wait(&neutrons_send_request,&status);

}

void send_fuel_assembly_neutrons(int count)
{
    int starting_index = fuel_assembly_neutrons_index[0];
    starting_index = 0;
    commit_sparse_neutrons_datatype(count); // Note: might have to change this for rolling messages

//    printf("Rank %d -- Sending neutrons...\n",rank);
    MPI_Issend(&neutrons[0],1,MPI_SPARSE_NEUTRONS, receiver, 0, world, &neutrons_send_request);
//    printf("Rank %d, count %ld\n",rank,count);
    MPI_Irecv(&neutrons[0],1,MPI_SPARSE_NEUTRONS,receiver,0,world,&neutrons_recv_request);
}







/**
 * Program entry point, this code is run with configuration file as command line argument
 **/
int main(int argc, char *argv[])
{

    MPI_Init(NULL,NULL);
    world = MPI_COMM_WORLD;
    MPI_Comm_rank(world,&rank);
    MPI_Comm_size(world,&size);

    sender = 1;
    receiver = 0;

    printf("Rank %d reporting\n",rank);

    if (argc != 3)
    {
        printf("You must provide two arguments, the configuration filename and output filename to the code\n");
        return -1;
    }

    time_t t;

    // Seed the random number generator
    srand((unsigned)time(&t));
    struct timeval start_time;


    // Parse the configuration and then initialise reactor core and neutrons from this
    struct simulation_configuration_struct configuration;
    parseConfiguration(argv[1], &configuration);
    initialiseReactorCore(&configuration);
    initialiseNeutrons(&configuration);

    commit_simple_neutron_datatype();
    commit_fission_event_datatype();

    // Empty the file we will use to store the reactor state
    clearReactorStateFile(argv[2]);
    if(rank==receiver) printf("Simulation configured for reactor core of size %dm by %dm by %dm, timesteps=%d dt=%dns\n", configuration.size_x,
                       configuration.size_y, configuration.size_z, configuration.num_timesteps, configuration.dt);
    if(rank==receiver) printf("------------------------------------------------------------------------------------------------\n");
    gettimeofday(&start_time, NULL); // Record simulation start time (for runtime statistics)

//    printf("Max int: %d. Num neutrons: %ld ",INT_MAX,configuration.max_neutrons);

    for (int i = 0; i < configuration.num_timesteps; i++)
    {
        // Progress in timesteps
        step(configuration.dt, &configuration);
        if (i > 0 && i % configuration.display_progess_frequency == 0)
        {
            if(rank==sender) calculateNumberActiveNeutrons(&configuration);
            if(rank==receiver) generateReport(configuration.dt, i, &configuration, start_time);
        }

        if (i > 0 && i % configuration.write_reactor_state_frequency == 0)
        {
            if(rank==receiver) writeReactorState(&configuration, i, argv[2]);
        }
    }


    // Now we are finished write some summary information
    unsigned long int num_fissions = getTotalNumberFissions(&configuration);
    double mev = getMeVFromFissions(num_fissions);
    double joules = getJoulesFromMeV(mev);
    if(rank==receiver) printf("------------------------------------------------------------------------------------------------\n");
    if(rank==receiver) printf("Model completed after %d timesteps\nTotal model time: %f secs\nTotal fissions: %ld releasing %e MeV and %e Joules\nTotal runtime: %.2f seconds\n",
                       configuration.num_timesteps, (NS_AS_SEC * configuration.dt) * configuration.num_timesteps, num_fissions, mev, joules, getElapsedTime(start_time));


    MPI_Finalize();

}

/**
 * Undertake a single timestep of processing
 **/
/*static*/ void step(int dt, struct simulation_configuration_struct *configuration)
{
    if(rank==sender) updateNeutrons(dt, configuration);
    if(rank==receiver) manageFuelAssemblyInteractions(configuration);
    if(rank==receiver) manageFuelAssemblyFissions(dt,configuration);
    if(rank==sender) updateReactorCore(dt, configuration);
}

/**
 * Writes a short report around the current state of the simulation to stdout
 **/
/*static*/ void generateReport(int dt, int timestep, struct simulation_configuration_struct *configuration, struct timeval start_time)
{
    unsigned long int num_active_neutrons = getNumberActiveNeutrons(configuration);
    unsigned long int num_fissions = getTotalNumberFissions(configuration);
    double mev = getMeVFromFissions(num_fissions);
    double joules = getJoulesFromMeV(mev);
    printf("Timestep: %d, model time is %e secs, current runtime is %.2f seconds. %ld active neutrons, %ld fissions, releasing %e MeV and %e Joules\n", timestep,
           (NS_AS_SEC * dt) * timestep, getElapsedTime(start_time), num_active_neutrons, num_fissions, mev, joules);
}

/**
 * Update the state of the reactor core at the current timestep, which will update the state
 * of each fuel assembly and neutron generator.
 **/
/*static*/ void updateReactorCore(int dt, struct simulation_configuration_struct *configuration)
{
    for (int i = 0; i < configuration->channels_x; i++)
    {
        for (int j = 0; j < configuration->channels_y; j++)
        {
            if (reactor_core[i][j].type == FUEL_ASSEMBLY)
            {
//                updateFuelAssembly(dt, &(reactor_core[i][j]));
                createNeutronsFromFission(&(reactor_core[i][j]));
            }

            if (reactor_core[i][j].type == NEUTRON_GENERATOR)
            {
                updateNeutronGenerator(dt, &(reactor_core[i][j]), configuration);
            }
        }
    }
}

/**
 * Update each neutron at the current timestep, moving its position based upon the velocity
 * components and energy, and then handling the collision of the neutron with a fuel channel,
 * moderator, or control rod
 **/
/*static*/ void updateNeutrons(int dt, struct simulation_configuration_struct *configuration)
{

    int num_fuel_interacting = 0; // Note: for rolling messages will probably have to move this

    for (long int i = 0; i < configuration->max_neutrons; i++)
    {
        if (neutrons[i].active)
        {
            updateNeutronPosition(&neutrons[i],dt);
            bool outbound = determineOutbound(&neutrons[i],configuration,i);

            struct channel_struct *reactorChannel = locateChannelFromPosition(neutrons[i].pos_x, neutrons[i].pos_y, configuration);

            if (reactorChannel != NULL && reactorChannel->type == FUEL_ASSEMBLY)
            {
                fuel_assembly_neutrons_index[num_fuel_interacting] = i;
                num_fuel_interacting++;
            }

        }
    }

    send_fuel_assembly_neutrons(num_fuel_interacting);

    for (long int i = 0; i < configuration->max_neutrons; i++)
    {
        if(neutrons[i].active)
        {
            // Now figure out if neutron is in a fuel assembly, moderator or control rod. If so then need to handle interaction
            struct channel_struct *reactorChannel = locateChannelFromPosition(neutrons[i].pos_x, neutrons[i].pos_y, configuration);
            interactWithReactor(&neutrons[i],reactorChannel,configuration,i);
        }
    }

    MPI_Status status;
    MPI_Wait(&neutrons_send_request,&status);
    MPI_Wait(&neutrons_recv_request,&status);
//    printf("Rank %d -- Received back.",rank);

    rebuildIndex(configuration);

}


/**
 * Update the state of a neutron generator for a timestep, generating the required
 * number of neutrons
 **/
/*static*/ void updateNeutronGenerator(int dt, struct channel_struct *channel, struct simulation_configuration_struct *configuration)
{
    unsigned long int number_new_neutrons = getNumberNeutronsFromGenerator(channel->contents.neutron_generator.weight, dt);
    for (int i = 0; i < number_new_neutrons; i++)
    {
        if (currentNeutronIndex == 0)
            break;
        currentNeutronIndex--;
        unsigned long int index = neutron_index[currentNeutronIndex];
        initialiseNeutron(&(neutrons[index]), channel, (double)(rand() / ((double)(RAND_MAX / configuration->size_z))));
    }
}

/**
 * Creates a specific number of neutrons that originate in a specific reactor channel
 * at a specific height (z) location
 **/
/*static*/ void createNeutrons(int num_neutrons, struct channel_struct *channel, double z)
{
    for (int k = 0; k < num_neutrons; k++)
    {
        if (currentNeutronIndex == 0)
            break;
        currentNeutronIndex--;
        unsigned long int index = neutron_index[currentNeutronIndex];
        initialiseNeutron(&(neutrons[index]), channel, z);
    }
}

/**
 * Initialises the reactor core at the start of the simulation from the configuration
 **/
/*static*/ void initialiseReactorCore(struct simulation_configuration_struct *simulation_configuration)
{

    max_fission_events = simulation_configuration->size_z / FUEL_PELLET_DEPTH;
    fission_array = (struct fission_event *) malloc(sizeof(struct fission_event) * max_fission_events); // Note: Comms related change

    reactor_core = (struct channel_struct **)malloc(sizeof(struct channel_struct *) * simulation_configuration->channels_x);

    for (int i = 0; i < simulation_configuration->channels_x; i++)
    {
        reactor_core[i] = (struct channel_struct *)malloc(sizeof(struct channel_struct) * simulation_configuration->channels_y);
    }
    for (int i = 0; i < simulation_configuration->channels_x; i++)
    {
        // Find the absolute x centre position of this channel
        double centre_x = ((i * CHANNEL_SIZE) + (CHANNEL_SIZE / 2));
        if (simulation_configuration->channel_layout_config == NULL)
        {
            // This channel has not been configured explicitly, hence assume it is empty
            for (int j = 0; j < simulation_configuration->channels_y; j++)
            {
                reactor_core[i][j].type = EMPTY_CHANNEL;
                reactor_core[i][j].x_centre = centre_x;
                reactor_core[i][j].y_centre = ((j * CHANNEL_SIZE) + (CHANNEL_SIZE / 2));
            }
        }
        else
        {
            for (int j = 0; j < simulation_configuration->num_channel_configs[i]; j++)
            {
                // For every configuration that was provided read what that was and initialise channel as required
                reactor_core[i][j].x_centre = centre_x;
                reactor_core[i][j].y_centre = ((j * CHANNEL_SIZE) + (CHANNEL_SIZE / 2));
                if (simulation_configuration->channel_layout_config[i][j] == CONFIG_MODERATOR)
                {
                    // This channel is a moderator, so set that and then the moderator type
                    reactor_core[i][j].type = MODERATOR;
                    if (simulation_configuration->moderator_type == WATER_MOD_TYPE_CONFIG)
                    {
                        reactor_core[i][j].contents.moderator.type = WATER;
                    }
                    else if (simulation_configuration->moderator_type == DEUTERIUM_MOD_TYPE_CONFIG)
                    {
                        reactor_core[i][j].contents.moderator.type = DEUTERIUM;
                    }
                    else if (simulation_configuration->moderator_type == GRAPHITE_MOD_TYPE_CONFIG)
                    {
                        reactor_core[i][j].contents.moderator.type = GRAPHITE;
                    }
                    else if (simulation_configuration->moderator_type == NONE_MOD_TYPE_CONFIG)
                    {
                        reactor_core[i][j].contents.moderator.type = NO_MODERATOR;
                    }
                }
                else if (simulation_configuration->channel_layout_config[i][j] == CONFIG_FUEL_ASSEMBLY)
                {
                    // This channel is a fuel assembly, so initialise that
                    reactor_core[i][j].type = FUEL_ASSEMBLY;
                    // Each fuel pellet is 40mm by 40mm by 2mm deep and weighs 1 gram
                    reactor_core[i][j].contents.fuel_assembly.num_pellets = simulation_configuration->size_z / FUEL_PELLET_DEPTH;
                    reactor_core[i][j].contents.fuel_assembly.num_fissions = 0;
                    reactor_core[i][j].contents.fuel_assembly.quantities = (double(*)[NUM_CHEMICALS])malloc(
                        sizeof(unsigned long int[NUM_CHEMICALS]) * reactor_core[i][j].contents.fuel_assembly.num_pellets);
                    for (int z = 0; z < reactor_core[i][j].contents.fuel_assembly.num_pellets; z++)
                    {
                        // For each pellet in the assembly set up the number of atoms present for each chemical, these
                        // will change as the simulation progresses and (hopefully!) fission occurs
                        for (int k = 0; k < NUM_CHEMICALS; k++)
                        {
                            enum chemical_type_enum chemical = getChemicalAtIndex(k);
                            if (chemical == UNKNOWN_CHEMICAL)
                            {
                                fprintf(stderr, "Unknown chemical at index '%d'\n", k);
                                exit(-1);
                            }
                            reactor_core[i][j].contents.fuel_assembly.quantities[z][k] = getAtomsPerGram(chemical) * (simulation_configuration->fuel_makeup_percentage[k] / 100.0);
                        }
                    }
                }
                else if (simulation_configuration->channel_layout_config[i][j] == CONFIG_CONTROL_ROD)
                {
                    // If the channel is a control rod then set this and store the absolute z location it is lowered to
                    reactor_core[i][j].type = CONTROL_ROD;
                    reactor_core[i][j].contents.control_rod.lowered_to_level = getControlRodLoweredToLevel(simulation_configuration, i, j);
                }
                else if (simulation_configuration->channel_layout_config[i][j] == CONFIG_EMPTY)
                {
                    reactor_core[i][j].type = EMPTY_CHANNEL;
                }
                else if (simulation_configuration->channel_layout_config[i][j] == CONFIG_NEUTRON_GENERATOR)
                {
                    reactor_core[i][j].type = NEUTRON_GENERATOR;
                    // Half a gram per cm in height
                    reactor_core[i][j].contents.neutron_generator.weight = simulation_configuration->size_z * 100 * NEUTRON_GENERATOR_WEIGHT_PER_CM;
                }
            }
            for (int j = simulation_configuration->num_channel_configs[i]; j < simulation_configuration->channels_y; j++)
            {
                // For any remaining channels that were not explicitly configured for this row, then set them as empty
                reactor_core[i][j].type = EMPTY_CHANNEL;
                reactor_core[i][j].x_centre = centre_x;
                reactor_core[i][j].y_centre = ((j * CHANNEL_SIZE) + (CHANNEL_SIZE / 2));
            }
        }
    }
}

/**
 * Initialises the neutron storage data structures at the start of the simulation so that all
 * neutrons are inactive. For performance we hold an array of indexes, each of which represents
 * empty slots in the neutron array which can be used. The currentNeutronIndex variable points
 * to the currnet end of the list, and when adding a neutron the value in that location of neutron_index
 * is read and currentNeutronIndex decremented. When deactivating a neutron the index of the newly freed
 * location is added to the next element of neutron_index and currentNeutronIndex is incremented
 **/
/*static*/ void initialiseNeutrons(struct simulation_configuration_struct *simulation_configuration)
{
    neutrons = (struct neutron_struct *)malloc(sizeof(struct neutron_struct) * simulation_configuration->max_neutrons);
    neutron_index = (unsigned long int *)malloc(sizeof(unsigned long int) * simulation_configuration->max_neutrons);
    fuel_assembly_neutrons_index = (int *) malloc(sizeof(int)*simulation_configuration->max_neutrons); // Note: the issue is not even message count, is indices for the indexed_block. Will need to improve datatype commit function to solve this
    for (int i = 0; i < simulation_configuration->max_neutrons; i++)
    {
        neutron_index[i] = i;
        fuel_assembly_neutrons_index[i] = -1; // Note: again, datatype issue could be problematic
        neutrons[i].active = false;
    }
    currentNeutronIndex = simulation_configuration->max_neutrons;
}

/**
 * For a control rod channel will return the absolute z height position that this has been lowered to based upon
 * the percentage insertion that was configured
 **/
/*static*/ double getControlRodLoweredToLevel(struct simulation_configuration_struct *simulation_configuration, int channel_x, int channel_y)
{
    int rodConfigurationIndex = findControlRodConfiguration(simulation_configuration, channel_x, channel_y);
    if (rodConfigurationIndex < 0)
    {
        fprintf(stderr, "Expected control rod configuration for x=%d y=%d but none can be found\n", channel_x, channel_y);
        exit(-1);
    }
    return simulation_configuration->size_z * (simulation_configuration->control_rod_configurations[rodConfigurationIndex].percentage / 100.0);
}

/**
 * Writes out the current state of the reactor at this timestep to a file
 **/
/*static*/ void writeReactorState(struct simulation_configuration_struct *configuration, int timestep, char *outfile)
{
    unsigned long int num_fissions = getTotalNumberFissions(configuration);
    double mev = getMeVFromFissions(num_fissions);
    double joules = getJoulesFromMeV(mev);
    FILE *f = fopen(outfile, "a");
    fprintf(f, "Reactor state at time %e secs, %ld fissions releasing %e MeV and %e Joules\n", (NS_AS_SEC * configuration->dt) * timestep, num_fissions, mev, joules);
    fprintf(f, "----------------------------------------------------------------------------\n");
    for (int i = 0; i < configuration->channels_x; i++)
    {
        for (int j = 0; j < configuration->channels_y; j++)
        {
            if (reactor_core[i][j].type == FUEL_ASSEMBLY)
            {
                double pc[11];
                getFuelAssemblyChemicalContents(&(reactor_core[i][j].contents.fuel_assembly), pc);
                fprintf(f, "Fuel assembly %d %d, %e U235 %e U238 %e Pu239 %e U236 %e Ba141 %e Kr92 %e Xe140 %e Sr94 %e Xe134 %e Zr103 %e Pu240\n",
                        i, j, pc[0], pc[1], pc[2], pc[3], pc[4], pc[5], pc[6], pc[7], pc[8], pc[9], pc[10]);
            }
        }
    }
    fprintf(f, "===========================================================================\n");
    fclose(f);
}

/**
 * Retrieves the quantities of atoms in a fuel assembly across all the pellets for
 * each chemical that will be written out to the file
 **/
/*static*/ void getFuelAssemblyChemicalContents(struct fuel_assembly_struct *fuel_rod, double *amounts)
{
    for (int i = 0; i < 11; i++)
        amounts[i] = 0;
    for (int i = 0; i < fuel_rod->num_pellets; i++)
    {
        for (int j = 0; j < 11; j++)
        {
            amounts[j] += fuel_rod->quantities[i][j];
        }
    }
}

/**
 * Clears out the file that we are going to write to for the reactor state, this is called at simulation
 * startup and it will overwrite any existing contents
 **/
/*static*/ void clearReactorStateFile(char *outfile)
{
    FILE *f = fopen(outfile, "w");
    fclose(f);
}

/**
 * Given an x and y position in the reactor core this will locate the channel
 * that that corresponds to
 **/
/*static*/ struct channel_struct *locateChannelFromPosition(double x, double y, struct simulation_configuration_struct *configuration)
{
    if (x > configuration->size_x || x < 0.0)
        return NULL;
    if (y > configuration->size_y || y < 0.0)
        return NULL;
    int channel_x = (int)(x / 0.2);
    int channel_y = (int)(y / 0.2);
    return &(reactor_core[channel_x][channel_y]);
}

/**
 * Based upon the properties of each fuel assembly will return the total number of fissions
 *  that have occured across all fuel assemblies in the simulation.
 **/
/*static*/ unsigned long int getTotalNumberFissions(struct simulation_configuration_struct *configuration)
{
    unsigned long int total_fissions = 0;
    for (int i = 0; i < configuration->channels_x; i++)
    {
        for (int j = 0; j < configuration->channels_y; j++)
        {
            if (reactor_core[i][j].type == FUEL_ASSEMBLY)
            {
                total_fissions += reactor_core[i][j].contents.fuel_assembly.num_fissions;
            }
        }
    }
    return total_fissions;
}


/**
 * Calculates the number of currently active neutrons in the simulation
 **/
/*static*/ unsigned long int calculateNumberActiveNeutrons(struct simulation_configuration_struct *configuration)
{

    unsigned long int activeNeutrons = 0;
    for (unsigned long int i = 0; i < configuration->max_neutrons; i++)
    {
        if (neutrons[i].active)
            activeNeutrons++;
    }

    MPI_Ssend(&activeNeutrons,1,MPI_UNSIGNED_LONG,receiver,1,world);

    return activeNeutrons;
}


/**
 * Determines the number of currently active neutrons in the simulation
 **/
/*static*/ unsigned long int getNumberActiveNeutrons(struct simulation_configuration_struct *configuration)
{

    MPI_Status status;
    unsigned long int activeNeutrons;

    MPI_Recv(&activeNeutrons,1,MPI_UNSIGNED_LONG,sender,1,world,&status); // Note: change to reduction using collective

    return activeNeutrons;
}


/**
 * Returns in seconds the elapsed time since the start_time argument and now
 **/
/*static*/ double getElapsedTime(struct timeval start_time)
{
    struct timeval curr_time;
    gettimeofday(&curr_time, NULL);
    long int elapsedtime = (curr_time.tv_sec * 1000000 + curr_time.tv_usec) - (start_time.tv_sec * 1000000 + start_time.tv_usec);
    return elapsedtime / 1000000.0;
}




void __attribute__ ((noinline)) updateNeutronPosition(struct neutron_struct *neutron, int dt) {

     // Rest mass is 1 for a neutron
    double total_velocity = MeVToVelocity(neutron->energy, 1);
    // These components are positive or negative which denote movement in one direction or another
    double component_velocity_x = (( neutron->x / 100.0) * total_velocity) * NS_AS_SEC * dt;
    double component_velocity_y = (( neutron->y / 100.0) * total_velocity) * NS_AS_SEC * dt;
    double component_velocity_z = (( neutron->z / 100.0) * total_velocity) * NS_AS_SEC * dt;

    neutron->pos_x += component_velocity_x;
    neutron->pos_y += component_velocity_y;
    neutron->pos_z += component_velocity_z;

}



/*
 * These functions were separated only for profiling ends. The compiler is inlining them all (which is desirable) but for testing they will be kept as noinline
 */

void __attribute__ ((noinline)) interactWithFuelAssembly(struct neutron_struct *neutron, struct channel_struct *reactorChannel, struct simulation_configuration_struct *configuration, long int i) {

    // It is in a fuel assembly channel, determine if it has collided with a neutron and if so deactivate it
    int fuel_pellet = (int)(neutron->pos_z / HEIGHT_FUEL_PELLET_M);
    if (fuel_pellet < reactorChannel->contents.fuel_assembly.num_pellets)
    {
        bool collision = determineAndHandleIfNeutronFuelCollision(neutron->energy, reactorChannel, fuel_pellet, configuration->collision_prob_multiplyer);
        if (collision)
        {
            neutron->active = false;
//            neutron_index[currentNeutronIndex] = i;
//            currentNeutronIndex++;
        }
    }
}


void __attribute__ ((noinline)) interactWithModerator(struct neutron_struct *neutron, struct channel_struct *reactorChannel, struct simulation_configuration_struct *configuration, long int i) {
    // The neutron is in the moderator, determine if it has been absorbed by the moderator or ot
    bool absorbed = determineAndHandleIfNeutronModeratorCollision(neutron, configuration->moderator_weight,
                                                                                  reactorChannel->contents.moderator.type, configuration->size_z);
    if (absorbed)
    {
        neutrons->active = false;
//        neutron_index[currentNeutronIndex] = i;
//        currentNeutronIndex++;
    }
}


void __attribute__ ((noinline)) interactWithControlRod(struct neutron_struct *neutron, struct channel_struct *reactorChannel, struct simulation_configuration_struct *configuration, long int i) {

    if (neutron->pos_z <= reactorChannel->contents.control_rod.lowered_to_level)
    {
        // Has hit the control rod, therefore this absorbed and removed from simulation
        neutron->active = false;
//        neutron_index[currentNeutronIndex] = i;
//        currentNeutronIndex++;
    }
}


void __attribute__ ((noinline)) interactWithReactor(struct neutron_struct *neutron, struct channel_struct *reactorChannel, struct simulation_configuration_struct *configuration, long int i) {

//    if (reactorChannel->type == FUEL_ASSEMBLY)
//    {
//        interactWithFuelAssembly(neutron,reactorChannel,configuration,i);
//    }

    if (reactorChannel->type == MODERATOR)
    {
        interactWithModerator(neutron,reactorChannel,configuration,i);
    }

    if (reactorChannel->type == CONTROL_ROD)
    {
        interactWithControlRod(neutron,reactorChannel,configuration,i);
    }

}

bool __attribute__ ((noinline)) determineOutbound(struct neutron_struct *neutron, struct simulation_configuration_struct *configuration, long int i) {

    bool outbound = neutron->pos_x > configuration->size_x || neutron->pos_x < 0.0 ||
           neutron->pos_y > configuration->size_y || neutron->pos_y < 0.0 ||
           neutron->pos_z > configuration->size_z || neutron->pos_z < 0.0 ;

    if(outbound) {
        // Moved out of the reactor core, so deactivate the neutron
        neutron->active = false;
//        neutron_index[currentNeutronIndex] = i;
//        currentNeutronIndex++;
    }

    return outbound;
}



void __attribute__ ((noinline)) rebuildIndex(struct simulation_configuration_struct *configuration)
{
    long int count_inactive = 0;
    currentNeutronIndex = configuration->max_neutrons;

    for (long int i = 0; i < configuration->max_neutrons; i++)
    {

        if (neutrons[i].active)
        {
            currentNeutronIndex --;
            neutron_index[currentNeutronIndex] = i;
        }
        else
        {
            neutron_index[count_inactive] = i;
            count_inactive ++;
        }

    }
}



/*
 * Parallelization-related functions
 */


void commit_sparse_neutrons_datatype(int count)
{
    MPI_Type_create_indexed_block(count,1,&fuel_assembly_neutrons_index[0],MPI_SIMPLE_NEUTRON,&MPI_SPARSE_NEUTRONS);
    MPI_Type_commit(&MPI_SPARSE_NEUTRONS);
}



void commit_simple_neutron_datatype()
{
    // Define block_lengths array.
    const int block_lengths[3] = {4,3,1};

    // Define displacements array.
    MPI_Aint displs[3];

    // Define dummy struct and get base address.
    struct neutron_struct dummy;
    MPI_Aint base_addr;
    MPI_Get_address(&dummy, &base_addr);

    //  Get addresses of dummy struct members and calculate displacements.
    MPI_Get_address(&dummy.pos_x, &displs[0]);
    MPI_Get_address(&dummy.x, &displs[1]);
    MPI_Get_address(&dummy.active, &displs[2]);
    displs[0] = MPI_Aint_diff(displs[0], base_addr);
    displs[1] = MPI_Aint_diff(displs[1], base_addr);
    displs[2] = MPI_Aint_diff(displs[2], base_addr);

    // Define types array.
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_SHORT, MPI_C_BOOL};

    // Create MPI datatype for neutron_struct using block_lengths, displs, and types arrays.
    MPI_Datatype tmp_type;
    MPI_Type_create_struct(3, block_lengths, displs, types, &tmp_type);
//    MPI_Type_create_struct(3, block_lengths, displs, types, &MPI_SIMPLE_NEUTRON);

    // Resize the MPI datatype to match the size of neutron_struct.
    MPI_Type_create_resized(tmp_type, 0, sizeof(dummy), &MPI_SIMPLE_NEUTRON);

    // Commit the MPI datatype.
    MPI_Type_commit(&MPI_SIMPLE_NEUTRON);

}

void commit_fission_event_datatype()
{
    // Define block_lengths array.
    const int block_lengths[2] = {1,1};

    // Define displacements array.
    MPI_Aint displs[2];

    // Define dummy struct and get base address.
    struct fission_event dummy;
    MPI_Aint base_addr;
    MPI_Get_address(&dummy, &base_addr);

    //  Get addresses of dummy struct members and calculate displacements.
    MPI_Get_address(&dummy.generated_neutrons, &displs[0]);
    MPI_Get_address(&dummy.pellet_height, &displs[1]);
    displs[0] = MPI_Aint_diff(displs[0], base_addr);
    displs[1] = MPI_Aint_diff(displs[1], base_addr);

    // Define types array.
    MPI_Datatype types[2] = {MPI_UNSIGNED_LONG, MPI_INT};

    // Create MPI datatype for neutron_struct using block_lengths, displs, and types arrays.
    MPI_Datatype tmp_type;
    MPI_Type_create_struct(2, block_lengths, displs, types, &tmp_type);

    // Resize the MPI datatype to match the size of neutron_struct.
    MPI_Type_create_resized(tmp_type, 0, sizeof(dummy), &MPI_FISSION_EVENT);

    // Commit the MPI datatype.
    MPI_Type_commit(&MPI_FISSION_EVENT);

}


//void commit_simple_neutron_datatype()
//{
//      // Define block_lengths array.
//      const int block_lengths[2] = {2,1};
//
//      // Define displacements array.
//      MPI_Aint displs[2];
//
//      // Define dummy struct and get base address.
//      struct neutron_struct dummy;
//      MPI_Aint base_addr;
//      MPI_Get_address(&dummy, &base_addr);
//
//      // Get addresses of dummy struct members and calculate displacements.
//      MPI_Get_address(&dummy.x, &displs[0]);
//      MPI_Get_address(&dummy.pos_x, &displs[1]);
//      displs[0] = MPI_Aint_diff(displs[0], base_addr);
//      displs[1] = MPI_Aint_diff(displs[1], base_addr);
//
//      // Define types array.
//      MPI_Datatype types[3] = {MPI_SHORT,MPI_DOUBLE,MPI_C_BOOL};
//
//      // Create MPI datatype for neutron_struct using block_lengths, displs, and types arrays.
//      MPI_Datatype tmp_type;
//      MPI_Type_create_struct(2, block_lengths, displs, types, &tmp_type);
////    MPI_Type_create_struct(2, block_lengths, displs, types, &MPI_SIMPLE_NEUTRON);
//
//      // Resize the MPI datatype to match the size of neutron_struct.
//      MPI_Type_create_resized(tmp_type, displs[0], sizeof(dummy), &MPI_SIMPLE_NEUTRON);
//
//      // Commit the MPI datatype.
//      MPI_Type_commit(&MPI_SIMPLE_NEUTRON);
//
//}



/**
 * Update the state of a specific fuel assembly in a channel for a timestep. This will fission all U236
 * and Pu239 in the assembly and update the constituent components as required
 **/
void executeFissions(int dt, struct channel_struct *channel)
{

    int initial_pellets = channel->contents.fuel_assembly.num_pellets;

    for (int i = 0; i < initial_pellets; i++)
    {

        fission_neutrons_counter = 0;

        unsigned long int num_u236 = (unsigned long int)channel->contents.fuel_assembly.quantities[i][U236];
        for (unsigned long int j = 0; j < num_u236; j++)
        {
            int num_neutrons = fissionU236(channel, i);
            createNeutrons(num_neutrons, channel, (i * HEIGHT_FUEL_PELLET_M) + (HEIGHT_FUEL_PELLET_M / 2));
            channel->contents.fuel_assembly.num_fissions++;
            fission_neutrons_counter+=num_neutrons;
        }
        unsigned long int num_pu240 = (unsigned long int)channel->contents.fuel_assembly.quantities[i][Pu240];
        for (unsigned long int j = 0; j < num_pu240; j++)
        {
            int num_neutrons = fissionPu240(channel, i);
            createNeutrons(num_neutrons, channel, (i * HEIGHT_FUEL_PELLET_M) + (HEIGHT_FUEL_PELLET_M / 2));
            channel->contents.fuel_assembly.num_fissions++;
            fission_neutrons_counter+=num_neutrons;
        }

        // Note: creates structure with i and the total number of neutrons and adds to an array.
        fission_array[i].generated_neutrons = fission_neutrons_counter;
        fission_array[i].pellet_height = i;

    }

    // Note: it is easiest to send (scatter) all fission_structs at once, once per channel (i.e. do in this function, it is required that every proc iterates channels in same order)
    MPI_Status status;
    MPI_Request request;
    MPI_Issend(&fission_array[0],initial_pellets,MPI_FISSION_EVENT,sender,0,world,&request);
    MPI_Wait(&request,&status);

}


void manageFuelAssemblyFissions(int dt, struct simulation_configuration_struct *configuration)
{
    for (int i = 0; i < configuration->channels_x; i++)
    {
        for (int j = 0; j < configuration->channels_y; j++)
        {
            if (reactor_core[i][j].type == FUEL_ASSEMBLY)
            {
                executeFissions(dt,&(reactor_core[i][j]));
            }
        }
    }
}




void createNeutronsFromFission(struct channel_struct *channel)
{
    MPI_Status status;
    int count;

    MPI_Recv(&fission_array[0],max_fission_events,MPI_FISSION_EVENT,receiver,0,world,&status); // Note: num_pellets on the neutron handler is not up-to-date, max_fission is needed
    MPI_Get_count(&status,MPI_FISSION_EVENT,&count);

    for(int i=0; i < count; i++)
    {
        createNeutrons(fission_array[i].generated_neutrons,channel,(fission_array[i].pellet_height * HEIGHT_FUEL_PELLET_M) + (HEIGHT_FUEL_PELLET_M / 2) );
    }

}
