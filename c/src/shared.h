// The neutrons that are currently moving throughout the reactor core
extern struct neutron_struct *neutrons;
// Indexes of empty neutrons (i.e. those inactive) that can be used
extern unsigned long int *neutron_index;
// The current index in the neutron_index array which is the next free neutron
extern unsigned long int currentNeutronIndex;
// The reactor core itself, each are channels in the x and y dimensions
extern struct channel_struct **reactor_core;



extern int size, rank;
extern MPI_Comm world;
extern MPI_Datatype MPI_SIMPLE_NEUTRON;
extern MPI_Datatype *MPI_SPARSE_NEUTRONS;
extern MPI_Datatype MPI_FISSION_EVENT;

extern MPI_Request *neutrons_send_request;
extern MPI_Request *neutrons_recv_request;

// Temporary variables

extern unsigned long int local_max_neutrons;
extern int *fuel_assembly_neutrons_index;
extern int fuelHandlerProc;
extern bool neutronHandler, fuelHandler,serial;
extern int max_fission_events;
extern int num_rolling_msgs;
extern int iterations_per_msg;
extern struct fission_event *fission_array;

