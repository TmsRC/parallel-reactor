#ifndef SERIAL_INCLUDE
#define SERIAL_INCLUDE
#include "simulation_configuration.h"
#include "simulation_support.h"

/*static*/ void serial_step(int, struct simulation_configuration_struct *);

/*static*/ void serial_updateReactorCore(int, struct simulation_configuration_struct *);
/*static*/ void serial_updateNeutrons(int, struct simulation_configuration_struct *);
/*static*/ void serial_updateFuelAssembly(int, struct channel_struct *);
/*static*/ void serial_updateNeutronGenerator(int, struct channel_struct *, struct simulation_configuration_struct *);
/*static*/ void serial_createNeutrons(int, struct channel_struct *, double);

/*static*/ unsigned long int serial_getNumberActiveNeutrons(struct simulation_configuration_struct *);

void serial_interactWithReactor(struct neutron_struct *neutron, struct channel_struct *reactorChannel, struct simulation_configuration_struct *configuration, long int i);

#endif
