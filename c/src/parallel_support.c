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
#include "parallel_support.h"


// The neutrons that are currently moving throughout the reactor core
struct neutron_struct *neutrons;
// Indexes of empty neutrons (i.e. those inactive) that can be used
unsigned long int *neutron_index;
// The current index in the neutron_index array which is the next free neutron
unsigned long int currentNeutronIndex = 0;
// The reactor core itself, each are channels in the x and y dimensions
struct channel_struct **reactor_core;

int global_counter=0;

// MPI variables
int size, rank;
MPI_Comm world;
MPI_Datatype MPI_SIMPLE_NEUTRON;
MPI_Datatype *MPI_SPARSE_NEUTRONS;
MPI_Datatype MPI_FISSION_EVENT;
MPI_Datatype MPI_STRIDED_FISSION;

MPI_Request *neutrons_send_request;
MPI_Request *neutrons_recv_request;

// Temporary variables
void commit_strided_fission_datatype(int);

unsigned long int local_max_neutrons;
// Note: I have to address this datatype issue
/*long*/ int *fuel_assembly_neutrons_index; // Contains the index in neutrons of neutrons that entered the fuel assembly. Also used to store the pellet height of generated neutrons

int fuelHandlerProc;
bool neutronHandler, fuelHandler;
int max_fission_events;
int num_rolling_msgs;
/*long*/ int iterations_per_msg;
struct fission_event *fission_array;



void manageFuelAssemblyInteractions(struct simulation_configuration_struct *configuration)
{

    MPI_Status status;

    int proc;
    int buff_loc;
    int proc_buff_size = local_max_neutrons/(size-1);

    int *received_messages = (int *) malloc(sizeof(int)*(size-1));
    int total_received_messages = 0;

    for(int i=0; i<(size-1); i++)
    {
        proc = i+1;
        buff_loc = i*proc_buff_size;
        MPI_Irecv(&neutrons[buff_loc],local_max_neutrons,MPI_SIMPLE_NEUTRON,proc,0,world,&neutrons_recv_request[i]); //Note: Need to change this so that I can use a simpler struct than neutrons in rank 0
        received_messages[i] = 0;
    }

    while( total_received_messages < num_rolling_msgs*(size-1) ) // Note: there should be exactly num_rolling * (size-1) messages to wait for
    {

        int i;
        MPI_Waitany(size-1,neutrons_recv_request,&i,&status);
        if(received_messages[i]>0) MPI_Wait(&neutrons_send_request[i],MPI_STATUS_IGNORE);

        received_messages[i] ++;
        total_received_messages ++;

        int count;
        MPI_Get_count(&status,MPI_SIMPLE_NEUTRON,&count);

        proc = i+1;
        buff_loc = i*proc_buff_size;


        for(long int j=0; j<count; j++)
        {
            struct channel_struct *reactorChannel = locateChannelFromPosition(neutrons[buff_loc+j].pos_x, neutrons[buff_loc+j].pos_y, configuration);
            interactWithFuelAssembly(&neutrons[buff_loc+j],reactorChannel,configuration,0);
        }

//        printf("Rank %d waiting to send to rank %d in manageFuelAssembly, message %d\n",rank,i,msg_num);
        MPI_Issend(&neutrons[buff_loc],count,MPI_SIMPLE_NEUTRON,proc,0,world,&neutrons_send_request[i]);
//        MPI_Wait(&neutrons_send_request[i],&status);


//        printf("Rank %d waiting to receive from rank %d in manageFuelAssembly, message %d\n",rank,i,msg_num);
        if(received_messages[i]<num_rolling_msgs) MPI_Irecv(&neutrons[buff_loc],local_max_neutrons,MPI_SIMPLE_NEUTRON,proc,0,world,&neutrons_recv_request[i]); //Note: Need to change this so that I can use a simpler struct than neutrons in rank 0
//        MPI_Wait(&neutrons_recv_request[i],&status);

    }

}


//void manageFuelAssemblyInteractions(struct simulation_configuration_struct *configuration)
//{
//
//    MPI_Status status;
//
//    int proc;
//    int buff_loc;
//    int proc_buff_size = local_max_neutrons/(size-1);
//
//    int *received_messages = (int *) malloc(sizeof(int)*(size-1));
//    int total_received_messages = 0;
//
//    for(int i=0; i<(size-1); i++)
//    {
//        proc = i+1;
//        buff_loc = i*proc_buff_size;
//        MPI_Irecv(&neutrons[buff_loc],local_max_neutrons,MPI_SIMPLE_NEUTRON,proc,0,world,&neutrons_recv_request[i]); //Note: Need to change this so that I can use a simpler struct than neutrons in rank 0
//        received_messages[i] = 0;
//    }
//
//    while( total_received_messages < num_rolling_msgs*(size-1) ) // Note: there should be exactly num_rolling * (size-1) messages to wait for
//    {
//        int received = 0;
//        for(int i=0; i<(size-1); i++)
//        {
//          MPI_Test(&neutrons_recv_request[i],&received,&status);
//          if(received)
//          {
//
//            received_messages[i] ++;
//            total_received_messages ++;
//
//            int count;
//            MPI_Get_count(&status,MPI_SIMPLE_NEUTRON,&count);
//
//            proc = i+1;
//            buff_loc = i*proc_buff_size;
//
//
//            for(long int j=0; j<count; j++)
//            {
//                struct channel_struct *reactorChannel = locateChannelFromPosition(neutrons[buff_loc+j].pos_x, neutrons[buff_loc+j].pos_y, configuration);
//                interactWithFuelAssembly(&neutrons[buff_loc+j],reactorChannel,configuration,0);
//            }
//
////            printf("Rank %d waiting to send to rank %d in manageFuelAssembly, message %d\n",rank,i,msg_num);
//            MPI_Issend(&neutrons[buff_loc],count,MPI_SIMPLE_NEUTRON,proc,0,world,&neutrons_send_request[i]);
////            MPI_Wait(&neutrons_send_request[i],&status);
//
//
////            printf("Rank %d waiting to receive from rank %d in manageFuelAssembly, message %d\n",rank,i,msg_num);
//            if(received_messages[i]<num_rolling_msgs) MPI_Irecv(&neutrons[buff_loc],local_max_neutrons,MPI_SIMPLE_NEUTRON,proc,0,world,&neutrons_recv_request[i]); //Note: Need to change this so that I can use a simpler struct than neutrons in rank 0
////            MPI_Wait(&neutrons_recv_request[i],&status);
//
//          }
//        }
//    }
//}


//void manageFuelAssemblyInteractions(struct simulation_configuration_struct *configuration)
//{
//
//    MPI_Status status;
//
//    int proc;
//    int buff_loc;
//    int proc_buff_size = local_max_neutrons/(size-1);
//
//
//    for(int msg_num = 0; msg_num < num_rolling_msgs; msg_num++) // Note: there should be exactly num_rolling * (size-1) messages to wait for
//    {
//        for(int i=0; i<(size-1); i++)
//        {
//
//            int count;
//            proc = i+1;
//            buff_loc = i*proc_buff_size;
////            printf("Rank %d waiting to receive from rank %d in manageFuelAssembly, message %d\n",rank,i,msg_num);
//            MPI_Irecv(&neutrons[buff_loc],local_max_neutrons,MPI_SIMPLE_NEUTRON,proc,0,world,&neutrons_recv_request[i]); //Note: Need to change this so that I can use a simpler struct than neutrons in rank 0
//            MPI_Wait(&neutrons_recv_request[i],&status);
//
//            MPI_Get_count(&status,MPI_SIMPLE_NEUTRON,&count);
//
//            for(long int j=0; j<count; j++)
//            {
//                struct channel_struct *reactorChannel = locateChannelFromPosition(neutrons[buff_loc+j].pos_x, neutrons[buff_loc+j].pos_y, configuration);
//                interactWithFuelAssembly(&neutrons[buff_loc+j],reactorChannel,configuration,0);
//            }
//
////            printf("Rank %d waiting to send to rank %d in manageFuelAssembly, message %d\n",rank,i,msg_num);
//            MPI_Issend(&neutrons[buff_loc],count,MPI_SIMPLE_NEUTRON,proc,0,world,&neutrons_send_request[i]);
//            MPI_Wait(&neutrons_send_request[i],&status);
//        }
//    }
//}

void createNeutronsFromFission(struct channel_struct *channel)
{
    MPI_Status status;
    int count;

    MPI_Recv(&fission_array[0],max_fission_events,MPI_FISSION_EVENT,fuelHandlerProc,0,world,&status); // Note: num_pellets on the neutron handler is not up-to-date, max_fission is needed
    MPI_Get_count(&status,MPI_FISSION_EVENT,&count);

    for(int i=0; i < count; i++)
    {
        createNeutrons(fission_array[i].generated_neutrons,channel,(fission_array[i].pellet_height * HEIGHT_FUEL_PELLET_M) + (HEIGHT_FUEL_PELLET_M / 2) );
    }

}



void executeFissions(int dt, struct channel_struct *channel)
{

    int initial_pellets = channel->contents.fuel_assembly.num_pellets;
    unsigned long int fission_neutrons_counter;

    for (int i = 0; i < initial_pellets; i++)
    {

        fission_neutrons_counter = 0;

        unsigned long int num_u236 = (unsigned long int)channel->contents.fuel_assembly.quantities[i][U236];
        for (unsigned long int j = 0; j < num_u236; j++)
        {
            int num_neutrons = fissionU236(channel, i);
            channel->contents.fuel_assembly.num_fissions++;
            fission_neutrons_counter+=num_neutrons;
        }
        unsigned long int num_pu240 = (unsigned long int)channel->contents.fuel_assembly.quantities[i][Pu240];
        for (unsigned long int j = 0; j < num_pu240; j++)
        {
            int num_neutrons = fissionPu240(channel, i);
            channel->contents.fuel_assembly.num_fissions++;
            fission_neutrons_counter+=num_neutrons;
        }

        // Note: creates structure with i and the total number of neutrons and adds to an array.
        fission_array[i].generated_neutrons = fission_neutrons_counter;
        fission_array[i].pellet_height = i;

    }

    sendGeneratedNeutrons(initial_pellets); // Note: this assumes global fission array

}


void sendGeneratedNeutrons(int initial_pellets)
{

    int num_receivers = size-1;

    int events_per_proc = initial_pellets/num_receivers;
    int remainder = initial_pellets%num_receivers;
    int displ;

    MPI_Status status;

    int start = 0;//rand(); // In order to distribute "fairly", but it seems to cause more problems than it solves

    commit_strided_fission_datatype(events_per_proc+1);
    for(int i=0; i<remainder; i++)
    {
        int proc = (start+i)%num_receivers+1;  // Note: assumes first proc is fuelHandler
        MPI_Issend(&fission_array[i],1,MPI_STRIDED_FISSION,proc,0,world,&neutrons_send_request[i]);
    }

    commit_strided_fission_datatype(events_per_proc);
    for(int i=remainder; i<size-1; i++)
    {
        int proc = (start+i)%num_receivers+1;
        MPI_Issend(&fission_array[i],1,MPI_STRIDED_FISSION,proc,0,world,&neutrons_send_request[i]);
    }



    for(int i=1; i<size; i++) // Note: this assumes fuelHandler = rank 0
    {
        MPI_Wait(&neutrons_send_request[i-1],&status);
    }

}





/**
 * Undertake a single timestep of processing
 **/
/*static*/ void step(int dt, struct simulation_configuration_struct *configuration)
{
    if(neutronHandler) updateNeutrons(dt, configuration);
    if(fuelHandler) manageFuelAssemblyInteractions(configuration);
    if(fuelHandler) manageFuelAssemblyFissions(dt,configuration);
    if(neutronHandler) updateReactorCore(dt, configuration);
}


/**
 * Update each neutron at the current timestep, moving its position based upon the velocity
 * components and energy, and then handling the collision of the neutron with a fuel channel,
 * moderator, or control rod
 **/
/*static*/ void updateNeutrons(int dt, struct simulation_configuration_struct *configuration)
{

    MPI_Status status;
    int num_fuel_interacting = 0; // Note: for rolling messages will probably have to move this
    int msg_num = 0;

    for (long int i = 0; i < local_max_neutrons; i++)
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

        // Note: lots of attention here, this loop is very different from the receiving one, so it has a significant deadlock risk 
        if ( ((i+1)%iterations_per_msg == 0 && (i+1)!=(num_rolling_msgs*iterations_per_msg)) || ((i+1)==local_max_neutrons) )
        {

            MPI_Status status;
            if(msg_num>0) MPI_Wait(&neutrons_send_request[msg_num-1],&status);
            if(msg_num>0) MPI_Wait(&neutrons_recv_request[msg_num-1],&status);

            sendFuelAssemblyNeutrons(num_fuel_interacting,msg_num);
            msg_num ++;
            num_fuel_interacting = 0;

        }

    }


    for (long int i = 0; i < local_max_neutrons; i++)
    {
        if(neutrons[i].active)
        {
            // Now figure out if neutron is in a fuel assembly, moderator or control rod. If so then need to handle interaction
            struct channel_struct *reactorChannel = locateChannelFromPosition(neutrons[i].pos_x, neutrons[i].pos_y, configuration);
            interactWithReactor(&neutrons[i],reactorChannel,configuration,i);
        }
    }

    MPI_Wait(&neutrons_send_request[msg_num-1],&status);
    MPI_Wait(&neutrons_recv_request[msg_num-1],&status);

    rebuildIndex(configuration);

}


/**
 * Update the state of a neutron generator for a timestep, generating the required
 * number of neutrons
 **/
/*static*/ void updateNeutronGenerator(int dt, struct channel_struct *channel, struct simulation_configuration_struct *configuration)
{
    unsigned long int number_new_neutrons = getNumberNeutronsFromGenerator(channel->contents.neutron_generator.weight, dt);
    number_new_neutrons = parallel_rescale_inv(number_new_neutrons,1,false,false);

    for (int i = 0; i < number_new_neutrons; i++)
    {
        if (currentNeutronIndex == 0)
            break;
        currentNeutronIndex--;
        unsigned long int index = neutron_index[currentNeutronIndex];
        initialiseNeutron(&(neutrons[index]), channel, (double)(rand() / ((double)(RAND_MAX / configuration->size_z))));
    }

}

void generateNeutrons(int dt, struct channel_struct *channel, struct simulation_configuration_struct *configuration)
{

    int num_remotes = size-1;
    unsigned long int *available_neutrons = (unsigned long int *) malloc(sizeof(unsigned long int)*size);
    int temp = 0;

    MPI_Gather(&temp,1,MPI_UNSIGNED_LONG,&available_neutrons,1,MPI_UNSIGNED_LONG,fuelHandlerProc,world);

    unsigned long int number_new_neutrons = getNumberNeutronsFromGenerator(channel->contents.neutron_generator.weight, dt);
    unsigned long int total_available = 0;
    unsigned long int total_allocated = 0;

    for(int i=0; i<size; i++)
    {
        total_available += available_neutrons[i-1];
    }

    for(int i=0; i<size && total_available!=0; i++)
    {
        available_neutrons[i-1] = (unsigned long int)( (double)(available_neutrons[i-1]) / total_available * number_new_neutrons);
        total_allocated = available_neutrons[i-1];
    }

    available_neutrons[size-1] += number_new_neutrons - total_allocated; // Note: assumes last proc is not fuelHandler

    MPI_Scatter(&available_neutrons,1,MPI_UNSIGNED_LONG,&temp,1,MPI_UNSIGNED_LONG,fuelHandlerProc,world);

    free(available_neutrons);

}


/**
 * Calculates the number of currently active neutrons in the simulation
 **/
/*static*/ unsigned long int calculateNumberActiveNeutrons(struct simulation_configuration_struct *configuration)
{
// It should be possible to calculate this as local_max_neutrons - currentNeutronIndex
    unsigned long int temp;
    unsigned long int activeNeutrons = 0;
    for (unsigned long int i = 0; i < local_max_neutrons; i++)
    {
        if (neutrons[i].active)
            activeNeutrons++;
    }

    MPI_Reduce(&activeNeutrons,&temp,1,MPI_UNSIGNED_LONG,MPI_SUM,fuelHandlerProc,world);

    return activeNeutrons;
}


/**
 * Determines the number of currently active neutrons in the simulation
 **/
/*static*/ unsigned long int getNumberActiveNeutrons(struct simulation_configuration_struct *configuration)
{
    unsigned long int activeNeutrons;
    unsigned long int temp = 0;

    MPI_Reduce(&temp,&activeNeutrons,1,MPI_UNSIGNED_LONG,MPI_SUM,fuelHandlerProc,world);

    return activeNeutrons;
}


void __attribute__ ((noinline)) interactWithReactor(struct neutron_struct *neutron, struct channel_struct *reactorChannel, struct simulation_configuration_struct *configuration, long int i) {

    if (reactorChannel->type == MODERATOR)
    {
        interactWithModerator(neutron,reactorChannel,configuration,i);
    }

    if (reactorChannel->type == CONTROL_ROD)
    {
        interactWithControlRod(neutron,reactorChannel,configuration,i);
    }

}

/*
 * Parallelization-related functions
 */


void commit_sparse_neutrons_datatype(int count, MPI_Datatype *sparse_neutrons)
{
    MPI_Type_create_indexed_block(count,1,&fuel_assembly_neutrons_index[0],MPI_SIMPLE_NEUTRON,sparse_neutrons);
    MPI_Type_commit(sparse_neutrons);
}



void commit_strided_fission_datatype(int count)
{
    MPI_Type_vector(count,1,(size-1),MPI_FISSION_EVENT,&MPI_STRIDED_FISSION);
    MPI_Type_commit(&MPI_STRIDED_FISSION);
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

//    commit_strided_fission_datatype();

}



// Ideally there would be multiple functions for different datatypes
unsigned long int parallel_rescale_inv(unsigned long int number, int num_ignored_processes, bool get_highest, bool get_total)
{
    if(get_total) return number;

    int involved_processes = size - num_ignored_processes;
    unsigned long int result = number/involved_processes;
    unsigned long int remainder = number%involved_processes;

    if(rank==(size-1) || get_highest) result += remainder;

    return result;

}




void sendFuelAssemblyNeutrons(int count, int msg_num)
{

    int starting_index = 0;
    commit_sparse_neutrons_datatype(count,&MPI_SPARSE_NEUTRONS[msg_num]); // Note: might have to change this for rolling messages


//    printf("Rank %d waiting to send to rank %d\n",rank,receiver);
    MPI_Issend(&neutrons[0],1,MPI_SPARSE_NEUTRONS[msg_num], fuelHandlerProc, 0, world, &neutrons_send_request[msg_num]);
    MPI_Irecv(&neutrons[0],1,MPI_SPARSE_NEUTRONS[msg_num],fuelHandlerProc,0,world,&neutrons_recv_request[msg_num]);
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
//                updateFuelAssembly(dt,&reactor_core[i][j]);
            }
        }
    }
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
                createNeutronsFromFission(&(reactor_core[i][j]));
//                createNeutronsFromFission_new(&(reactor_core[i][j]));
            }

            if (reactor_core[i][j].type == NEUTRON_GENERATOR)
            {
                updateNeutronGenerator(dt, &(reactor_core[i][j]), configuration);
            }
        }
    }
}
