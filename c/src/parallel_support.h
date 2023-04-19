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
void commit_sparse_neutrons_datatype(int count, MPI_Datatype *sparse_neutrons);
unsigned long int parallel_rescale_inv(unsigned long int number, int num_ignored_processes, bool get_highest, bool get_total);

void executeFissions(int dt, struct channel_struct *channel);
void manageFuelAssemblyFissions(int dt, struct simulation_configuration_struct *configuration);
void createNeutronsFromFission(struct channel_struct *channel);
unsigned long int calculateNumberActiveNeutrons(struct simulation_configuration_struct *);

void sendFuelAssemblyNeutrons(int count, int msg_num);
void sendGeneratedNeutrons(int initial_pellets);
