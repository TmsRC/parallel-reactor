#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "simulation_configuration.h"
#include "simulation_support.h"
#include "serial_support.h"
#include "parallel_support.h"
#include "shared.h"

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
//
//// The neutrons that are currently moving throughout the reactor core
//struct neutron_struct *neutrons;
//// Indexes of empty neutrons (i.e. those inactive) that can be used
//unsigned long int *neutron_index;
//// The current index in the neutron_index array which is the next free neutron
//unsigned long int currentNeutronIndex = 0;
//// The reactor core itself, each are channels in the x and y dimensions
//struct channel_struct **reactor_core;
//
///*static*/ void serial_step(int, struct simulation_configuration_struct *);
//
///*static*/ void generateReport(int, int, struct simulation_configuration_struct *, struct timeval);
//
///*static*/ void serial_updateReactorCore(int, struct simulation_configuration_struct *);
///*static*/ void serial_updateNeutrons(int, struct simulation_configuration_struct *);
///*static*/ void serial_updateFuelAssembly(int, struct channel_struct *);
///*static*/ void serial_updateNeutronGenerator(int, struct channel_struct *, struct simulation_configuration_struct *);
///*static*/ void serial_createNeutrons(int, struct channel_struct *, double);
//
///*static*/ void initialiseReactorCore(struct simulation_configuration_struct *);
///*static*/ void initialiseNeutrons(struct simulation_configuration_struct *);
///*static*/ double getControlRodLoweredToLevel(struct simulation_configuration_struct *, int, int);
///*static*/ void writeReactorState(struct simulation_configuration_struct *, int, char *);
///*static*/ void getFuelAssemblyChemicalContents(struct fuel_assembly_struct *, double *);
///*static*/ void clearReactorStateFile(char *);
///*static*/ struct channel_struct *locateChannelFromPosition(double, double, struct simulation_configuration_struct *);
///*static*/ unsigned long int getTotalNumberFissions(struct simulation_configuration_struct *);
//
///*static*/ unsigned long int serial_getNumberActiveNeutrons(struct simulation_configuration_struct *);
//
///*static*/ double getElapsedTime(struct timeval);
//
//void updateNeutronPosition(struct neutron_struct *neutron, int dt);
//void interactWithFuelAssembly(struct neutron_struct *neutron, struct channel_struct *reactorChannel, struct simulation_configuration_struct *configuration, long int i);
//void interactWithModerator(struct neutron_struct *neutron, struct channel_struct *reactorChannel, struct simulation_configuration_struct *configuration, long int i);
//void interactWithControlRod(struct neutron_struct *neutron, struct channel_struct *reactorChannel, struct simulation_configuration_struct *configuration, long int i);
//
//void serial_interactWithReactor(struct neutron_struct *neutron, struct channel_struct *reactorChannel, struct simulation_configuration_struct *configuration, long int i);
//
//bool determineOutbound(struct neutron_struct *neutron, struct simulation_configuration_struct *configuration, long int i);


/**
 * Undertake a single timestep of processing
 **/
/*static*/ void serial_step(int dt, struct simulation_configuration_struct *configuration)
{
    serial_updateNeutrons(dt, configuration);
    serial_updateReactorCore(dt, configuration);
}

/**
 * Update the state of the reactor core at the current timestep, which will update the state
 * of each fuel assembly and neutron generator.
 **/
/*static*/ void serial_updateReactorCore(int dt, struct simulation_configuration_struct *configuration)
{
    for (int i = 0; i < configuration->channels_x; i++)
    {
        for (int j = 0; j < configuration->channels_y; j++)
        {
            if (reactor_core[i][j].type == FUEL_ASSEMBLY)
            {
                serial_updateFuelAssembly(dt, &(reactor_core[i][j]));
            }

            if (reactor_core[i][j].type == NEUTRON_GENERATOR)
            {
                serial_updateNeutronGenerator(dt, &(reactor_core[i][j]), configuration);
            }
        }
    }
}

/**
 * Update each neutron at the current timestep, moving its position based upon the velocity
 * components and energy, and then handling the collision of the neutron with a fuel channel,
 * moderator, or control rod
 **/
/*static*/ void serial_updateNeutrons(int dt, struct simulation_configuration_struct *configuration)
{
    for (long int i = 0; i < configuration->max_neutrons; i++)
    {
        if (neutrons[i].active)
        {
            updateNeutronPosition(&neutrons[i],dt);
            bool outbound = determineOutbound(&neutrons[i],configuration,i);
        }
    }

    for (long int i = 0; i < configuration->max_neutrons; i++)
    {
        if(neutrons[i].active)
        {

            // Now figure out if neutron is in a fuel assembly, moderator or control rod. If so then need to handle interaction
            struct channel_struct *reactorChannel = locateChannelFromPosition(neutrons[i].pos_x, neutrons[i].pos_y, configuration);
            if (reactorChannel != NULL)
            {
                serial_interactWithReactor(&neutrons[i],reactorChannel,configuration,i);
            }
            else
            {
                fprintf(stderr, "Unable to locate reactor core channel for x=%f and y=%f\n", neutrons[i].pos_x, neutrons[i].pos_y);
                exit(-1);
            }
        }
    }

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

/**
 * Update the state of a specific fuel assembly in a channel for a timestep. This will fission all U236
 * and Pu239 in the assembly and update the constituent components as required
 **/
/*static*/ void serial_updateFuelAssembly(int dt, struct channel_struct *channel)
{
    for (int i = 0; i < channel->contents.fuel_assembly.num_pellets; i++)
    {
        unsigned long int num_u236 = (unsigned long int)channel->contents.fuel_assembly.quantities[i][U236];
        for (unsigned long int j = 0; j < num_u236; j++)
        {
            int num_neutrons = fissionU236(channel, i);
            createNeutrons(num_neutrons, channel, (i * HEIGHT_FUEL_PELLET_M) + (HEIGHT_FUEL_PELLET_M / 2));
            channel->contents.fuel_assembly.num_fissions++;
        }
        unsigned long int num_pu240 = (unsigned long int)channel->contents.fuel_assembly.quantities[i][Pu240];
        for (unsigned long int j = 0; j < num_pu240; j++)
        {
            int num_neutrons = fissionPu240(channel, i);
            createNeutrons(num_neutrons, channel, (i * HEIGHT_FUEL_PELLET_M) + (HEIGHT_FUEL_PELLET_M / 2));
            channel->contents.fuel_assembly.num_fissions++;
        }
    }
}

/**
 * Update the state of a neutron generator for a timestep, generating the required
 * number of neutrons
 **/
/*static*/ void serial_updateNeutronGenerator(int dt, struct channel_struct *channel, struct simulation_configuration_struct *configuration)
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
 * Determines the number of currently active neutrons in the simulation
 **/
/*static*/ unsigned long int serial_getNumberActiveNeutrons(struct simulation_configuration_struct *configuration)
{
    unsigned long int activeNeutrons = 0;
    for (unsigned long int i = 0; i < configuration->max_neutrons; i++)
    {
        if (neutrons[i].active)
            activeNeutrons++;
    }
    return activeNeutrons;
}


void __attribute__ ((noinline)) serial_interactWithReactor(struct neutron_struct *neutron, struct channel_struct *reactorChannel, struct simulation_configuration_struct *configuration, long int i) {

    if (reactorChannel->type == FUEL_ASSEMBLY)
    {
        interactWithFuelAssembly(neutron,reactorChannel,configuration,i);
    }

    if (reactorChannel->type == MODERATOR)
    {
        interactWithModerator(neutron,reactorChannel,configuration,i);
    }

    if (reactorChannel->type == CONTROL_ROD)
    {
        interactWithControlRod(neutron,reactorChannel,configuration,i);
    }

}
