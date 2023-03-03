module reactor_simulator
  use configuration, only : simulation_configuration_type, parseConfiguration, findControlRodConfiguration, &
      CONFIG_EMPTY, CONFIG_FUEL_ASSEMBLY, CONFIG_CONTROL_ROD, CONFIG_NEUTRON_GENERATOR, CONFIG_MODERATOR, &
      WATER_MOD_TYPE_CONFIG, DEUTERIUM_MOD_TYPE_CONFIG, GRAPHITE_MOD_TYPE_CONFIG, NONE_MOD_TYPE_CONFIG
  use simulation_support, only : channel_type, neutron_type, fuel_assembly_type, &
                                getAtomsPerGram, MeVToVelocity, fissionU236, fissionPu240, initialiseNeutron, &
                                determineAndHandleIfNeutronFuelCollision, determineAndHandleIfNeutronModeratorCollision, &
                                getNumberNeutronsFromGenerator, getMeVFromFissions, getJoulesFromMeV, &
                                seedRandomNumberGenerator, getRandomReal, getRandomInt, Pu240, U236, &
                                WATER, DEUTERIUM, GRAPHITE, NO_MODERATOR, MODERATOR, FUEL_ASSEMBLY, CONTROL_ROD, &
                                EMPTY_CHANNEL, NEUTRON_GENERATOR, NUM_CHEMICALS, UNKNOWN_CHEMICAL

  implicit none

  ! Size of each reactor channel in meters (they are cuboid, so this value in x and y)
  real(kind=8), parameter :: CHANNEL_SIZE=0.2
  ! Height of a fuel pellet in meters (they are 40mm by 40mm by 2.5mm)
  real(kind=8), parameter :: HEIGHT_FUEL_PELLET_M=0.0025
  ! Depth of fuel pellets in meters (they are x=40mm by y=40mm by z=2mm)
  real(kind=8), parameter :: FUEL_PELLET_DEPTH=0.002
  ! Weight in grams of the neutron generator for each cm in height
  real(kind=8), parameter :: NEUTRON_GENERATOR_WEIGHT_PER_CM=0.5
  ! Number of nanoseconds in a second
  real(kind=8), parameter :: NS_AS_SEC=1e-9

  ! The reactor core itself, each are channels in the x and y dimensions
  type(channel_type), dimension(:,:), allocatable :: reactor_core
  type(neutron_type), dimension(:), allocatable :: neutrons
  integer*8, dimension(:), allocatable :: neutron_index
  integer*8 :: currentNeutronIndex
contains

  ! Program entry point, this code is run with configuration file as command line argument
  subroutine run_simulation(config_file, output_file)
    character(len=*), intent(in) :: config_file, output_file

    type(simulation_configuration_type) :: sim_config
    integer :: start_time_values(8), i, end_time_values(8)
    integer*8 :: num_fissions
    real(kind=8) :: mev, joules

    call seedRandomNumberGenerator()
    ! Parse the configuration and then initialise reactor core and neutrons from this
    call parseConfiguration(config_file, sim_config)
    call initialiseReactorCore(sim_config)
    call initialiseNeutrons(sim_config)
    ! Empty the file we will use to store the reactor state
    call clearReactorStateFile(output_file)
    print *, "Simulation configured for reactor core of size ", &
      trim(adjustl(intToString(sim_config%size_x))), &
      " by ", trim(adjustl(intToString(sim_config%size_y))), " by ", &
      trim(adjustl(intToString(sim_config%size_z))), ", timesteps=", &
      trim(adjustl(intToString(sim_config%num_timesteps))), " dt=", &
      trim(adjustl(intToString(sim_config%dt)))
    ! Record simulation start time (for runtime statistics)
    call date_and_time(values=start_time_values)
    do i=1, sim_config%num_timesteps
      ! Progress in timesteps
      call step(sim_config%dt, sim_config)

      if (i .gt. 0 .and. mod(i, sim_config%display_progess_frequency) == 0) then
        call generateReport(sim_config%dt, i, sim_config, start_time_values)
      end if

      if (i .gt. 0 .and. mod(i, sim_config%write_reactor_state_frequency) == 0) then
        call writeReactorState(sim_config, i, output_file)
      end if
    end do

    ! Now we are finished write some summary information
    num_fissions=getTotalNumberFissions(sim_config)
    mev=getMeVFromFissions(num_fissions)
    joules=getJoulesFromMeV(mev)
    call date_and_time(values=end_time_values)
    print *, "------------------------------------------------------------------------------------------------"
    print *, "Model completed after ", trim(adjustl(intToString(sim_config%num_timesteps))), " timesteps"
    print *, "Total model time: ", trim(adjustl(realToExprString((NS_AS_SEC*sim_config%dt)*sim_config%num_timesteps)))
    print *, "Total fissions: ", trim(adjustl(longToString(num_fissions))), " fissions releasing ", &
        trim(adjustl(realToExprString(mev)))," MeV and ", trim(adjustl(realToExprString(joules))), " Joules"
    print *, "Total runtime: ", trim(adjustl(realToString(getTimeDifference(start_time_values, end_time_values)))), " seconds"
  end subroutine run_simulation

  ! Undertake a single timestep of processing
  subroutine step(dt, sim_config)
    integer, intent(in) :: dt
    type(simulation_configuration_type), intent(in) :: sim_config

    call updateNeutrons(dt, sim_config)
    call updateReactorCore(dt, sim_config)
  end subroutine step

  ! Writes a short report around the current state of the simulation to stdout
  subroutine generateReport(dt, timestep, sim_config, starttime)
    integer, intent(in) :: dt, timestep, starttime(8)
    type(simulation_configuration_type), intent(in) :: sim_config

    integer*8 :: num_fissions
    real(kind=8) :: mev, joules
    integer :: end_time_values(8)

    num_fissions=getTotalNumberFissions(sim_config)
    mev=getMeVFromFissions(num_fissions)
    joules=getJoulesFromMeV(mev)

    call date_and_time(values=end_time_values)

    print *, "Timestep: ", trim(adjustl(intToString(timestep))), ", model time is ", &
      trim(adjustl(realToExprString((NS_AS_SEC*dt)*timestep))), " sec, current runtime is ", &
      trim(adjustl(realToString(getTimeDifference(starttime, end_time_values)))), " seconds. ", &
      trim(adjustl(longToString(getNumberActiveNeutrons(sim_config)))), " active neutrons, ", &
      trim(adjustl(longToString(num_fissions))), " fissions, releasing ", &
      trim(adjustl(realToExprString(mev)))," MeV and ", trim(adjustl(realToExprString(joules))), " Joules"
  end subroutine generateReport

  ! Update the state of the reactor core at the current timestep, which will update the state
  ! of each fuel assembly and neutron generator.
  subroutine updateReactorCore(dt, sim_config)
    integer, intent(in) :: dt
    type(simulation_configuration_type), intent(in) :: sim_config

    integer :: i, j

    do i=0, sim_config%channels_x-1
      do j=0, sim_config%channels_y-1
        if (reactor_core(i,j)%type_of_channel == FUEL_ASSEMBLY) then
          call updateFuelAssembly(dt, reactor_core(i,j))
        else if (reactor_core(i,j)%type_of_channel == NEUTRON_GENERATOR) then
          call updateNeutronGenerator(dt, reactor_core(i,j), sim_config)
        end if
      end do
    end do
  end subroutine updateReactorCore

  ! Update each neutron at the current timestep, moving its position based upon the velocity
  ! components and energy, and then handling the collision of the neutron with a fuel channel,
  ! moderator, or control rod
  subroutine updateNeutrons(dt, sim_config)
    integer, intent(in) :: dt
    type(simulation_configuration_type), intent(in) :: sim_config

    integer i, channel_x, channel_y, fuel_pellet
    logical :: collision, absorbed
    real(kind=8) :: total_velocity, component_velocity_x, component_velocity_y, component_velocity_z

    do i=1, sim_config%max_neutrons
      if (neutrons(i)%active) then
        ! Rest mass is 1 for a neutron
        total_velocity=MeVToVelocity(neutrons(i)%energy, 1)
        ! These components are positive or negative which denote movement in one direction or another
        component_velocity_x=((abs(neutrons(i)%x)/100.0) * total_velocity)*NS_AS_SEC*dt
        component_velocity_y=((abs(neutrons(i)%y)/100.0) * total_velocity)*NS_AS_SEC*dt
        component_velocity_z=((abs(neutrons(i)%z)/100.0) * total_velocity)*NS_AS_SEC*dt
        if (neutrons(i)%x .gt. 0) then
          neutrons(i)%pos_x=neutrons(i)%pos_x+component_velocity_x
        else
          neutrons(i)%pos_x=neutrons(i)%pos_x-component_velocity_x
        end if
        if (neutrons(i)%y .gt. 0) then
          neutrons(i)%pos_y=neutrons(i)%pos_y+component_velocity_y
        else
          neutrons(i)%pos_y=neutrons(i)%pos_y-component_velocity_y
        end if
        if (neutrons(i)%z .gt. 0) then
          neutrons(i)%pos_z=neutrons(i)%pos_z+component_velocity_z
        else
          neutrons(i)%pos_z=neutrons(i)%pos_z-component_velocity_z
        end if

        if (neutrons(i)%pos_x .gt. sim_config%size_x .or. neutrons(i)%pos_x .lt. 0.0 &
            .or. neutrons(i)%pos_y .gt. sim_config%size_y .or. neutrons(i)%pos_y .lt. 0.0 &
            .or. neutrons(i)%pos_z .gt. sim_config%size_z .or. neutrons(i)%pos_z .lt. 0.0) then
          ! Moved out of the reactor core, so deactivate the neutron
          neutrons(i)%active=.false.
          neutron_index(currentNeutronIndex)=i
          currentNeutronIndex=currentNeutronIndex+1
          cycle
        end if
        ! Now figure out if neutron is in a fuel assembly, moderator or control rod. If so then need to handle interaction
        if (locateChannelFromPosition(neutrons(i)%pos_x, neutrons(i)%pos_y, sim_config, channel_x, channel_y)) then
          if (reactor_core(channel_x, channel_y)%type_of_channel == FUEL_ASSEMBLY) then
            ! It is in a fuel assembly channel, determine if it has collided with a neutron and if so deactivate it
            fuel_pellet=neutrons(i)%pos_z / HEIGHT_FUEL_PELLET_M
            if (fuel_pellet .lt. reactor_core(channel_x, channel_y)%fuel_assembly%num_pellets) then
              collision=determineAndHandleIfNeutronFuelCollision(neutrons(i)%energy, reactor_core(channel_x, channel_y), &
                fuel_pellet, sim_config%collision_prob_multiplyer)
              if (collision) then
                neutrons(i)%active=.false.
                neutron_index(currentNeutronIndex)=i
                currentNeutronIndex=currentNeutronIndex+1
              end if
            end if
          else if (reactor_core(channel_x, channel_y)%type_of_channel == MODERATOR) then
            ! The neutron is in the moderator, determine if it has been absorbed by the moderator or ot
            absorbed=determineAndHandleIfNeutronModeratorCollision(neutrons(i), sim_config%moderator_weight, &
                              reactor_core(channel_x, channel_y)%moderator%moderator_material, sim_config%size_z)
            if (absorbed) then
              neutrons(i)%active=.false.
              neutron_index(currentNeutronIndex)=i
              currentNeutronIndex=currentNeutronIndex+1
            end if
          else if (reactor_core(channel_x, channel_y)%type_of_channel == CONTROL_ROD) then
            if (neutrons(i)%pos_z .le. reactor_core(channel_x, channel_y)%control_rod%lowered_to_level) then
              ! Has hit the control rod, therefore this absorbed and removed from simulation
              neutrons(i)%active=.false.
              neutron_index(currentNeutronIndex)=i
              currentNeutronIndex=currentNeutronIndex+1
            end if
          end if
        else
          print *, "Unable to locate reactor core channel for x=", neutrons(i)%pos_x, " and y=", neutrons(i)%pos_y
          call exit(-1)
        end if
      end if
    end do
  end subroutine updateNeutrons

  ! Update the state of a specific fuel assembly in a channel for a timestep. This will fission all U236
  ! and Pu239 in the assembly and update the constituent components as required
  subroutine updateFuelAssembly(dt, reactor_channel)
    integer, intent(in) :: dt
    type(channel_type), intent(inout) :: reactor_channel

    integer :: i, num_neutrons
    integer*8 :: j, num_u236, num_pu240

    do i=1, reactor_channel%fuel_assembly%num_pellets
      num_u236=int(reactor_channel%fuel_assembly%quantities(U236, i), 8)
      do j=1, num_u236
        num_neutrons=fissionU236(reactor_channel, i)
        call createNeutrons(num_neutrons, reactor_channel, (i*HEIGHT_FUEL_PELLET_M)+(HEIGHT_FUEL_PELLET_M/2))
        reactor_channel%fuel_assembly%num_fissions=&
          reactor_channel%fuel_assembly%num_fissions+1
      end do
      num_pu240=int(reactor_channel%fuel_assembly%quantities(Pu240, i), 8)
      do j=1, num_pu240
        num_neutrons=fissionPu240(reactor_channel, i)
        call createNeutrons(num_neutrons, reactor_channel, (i*HEIGHT_FUEL_PELLET_M)+(HEIGHT_FUEL_PELLET_M/2))
        reactor_channel%fuel_assembly%num_fissions=&
          reactor_channel%fuel_assembly%num_fissions+1
      end do
    end do
  end subroutine updateFuelAssembly

  ! Update the state of a neutron generator for a timestep, generating the required
  ! number of neutrons
  subroutine updateNeutronGenerator(dt, reactor_channel, sim_config)
    integer, intent(in) :: dt
    type(channel_type), intent(inout) :: reactor_channel
    type(simulation_configuration_type), intent(in) :: sim_config

    integer :: i
    integer *8 :: number_new_neutrons, idx

    number_new_neutrons=getNumberNeutronsFromGenerator(reactor_channel%neutron_generator%weight, dt)
    do i=1, number_new_neutrons
      if (currentNeutronIndex == 1) exit
      currentNeutronIndex=currentNeutronIndex-1
      idx=neutron_index(currentNeutronIndex)
      call initialiseNeutron(neutrons(idx), reactor_channel, real(getRandomInt(sim_config%size_z*1000) / 1000.0, 8))
    end do
  end subroutine updateNeutronGenerator

  ! Creates a specific number of neutrons that originate in a specific reactor channel
  ! at a specific height (z) location
  subroutine createNeutrons(num_neutrons, reactor_channel, z)
    integer, intent(in) :: num_neutrons
    type(channel_type), intent(inout) :: reactor_channel
    real(kind=8), intent(in) :: z

    integer :: k

    do k=1, num_neutrons
      if (currentNeutronIndex == 1) exit
      currentNeutronIndex=currentNeutronIndex-1
      call initialiseNeutron(neutrons(neutron_index(currentNeutronIndex)), reactor_channel, z)
    end do
  end subroutine createNeutrons

  ! Initialises the reactor core at the start of the simulation from the configuration
  subroutine initialiseReactorCore(sim_config)
    type(simulation_configuration_type), intent(in) :: sim_config

    integer :: i, j, k, z
    real(kind=8) :: centre_x

    allocate(reactor_core(0:sim_config%channels_x-1, 0:sim_config%channels_y-1))
    do i=0, sim_config%channels_x-1
      !  Find the absolute x centre position of this channel
      centre_x=((i*CHANNEL_SIZE) + (CHANNEL_SIZE/2))
      do j=0, sim_config%num_channel_configs(i)-1
        ! For every configuration that was provided read what that was and initialise channel as required
        reactor_core(i,j)%x_centre=centre_x
        reactor_core(i,j)%y_centre=((j*CHANNEL_SIZE) + (CHANNEL_SIZE/2))
        if (sim_config%channel_layout_config(i,j) == CONFIG_MODERATOR) then
          ! This channel is a moderator, so set that and then the moderator type
          reactor_core(i,j)%type_of_channel=MODERATOR
          if (sim_config%moderator_type == WATER_MOD_TYPE_CONFIG) then
            reactor_core(i,j)%moderator%moderator_material=WATER
          else if (sim_config%moderator_type == DEUTERIUM_MOD_TYPE_CONFIG) then
            reactor_core(i,j)%moderator%moderator_material=DEUTERIUM
          else if (sim_config%moderator_type == GRAPHITE_MOD_TYPE_CONFIG) then
            reactor_core(i,j)%moderator%moderator_material=GRAPHITE
          else if (sim_config%moderator_type == NONE_MOD_TYPE_CONFIG) then
            reactor_core(i,j)%moderator%moderator_material=NO_MODERATOR
          end if
        else if (sim_config%channel_layout_config(i,j) == CONFIG_FUEL_ASSEMBLY) then
          ! This channel is a fuel assembly, so initialise that
          reactor_core(i,j)%type_of_channel=FUEL_ASSEMBLY
          ! Each fuel pellet is 40mm by 40mm by 2mm deep and weighs 1 gram
          reactor_core(i,j)%fuel_assembly%num_pellets=sim_config%size_z / FUEL_PELLET_DEPTH
          reactor_core(i,j)%fuel_assembly%num_fissions=0
          allocate(reactor_core(i,j)%fuel_assembly%quantities(NUM_CHEMICALS, &
            reactor_core(i,j)%fuel_assembly%num_pellets))
          do k=1, reactor_core(i,j)%fuel_assembly%num_pellets
            ! For each pellet in the assembly set up the number of atoms present for each chemical, these
            ! will change as the simulation progresses and (hopefully!) fission occurs
            do z=1, NUM_CHEMICALS
              reactor_core(i,j)%fuel_assembly%quantities(z, k)=getAtomsPerGram(z) * (sim_config%fuel_makeup_percentage(z) / 100.0)
            end do
          end do
        else if (sim_config%channel_layout_config(i,j) == CONFIG_CONTROL_ROD) then
          ! If the channel is a control rod then set this and store the absolute z location it is lowered to
          reactor_core(i,j)%type_of_channel=CONTROL_ROD
          reactor_core(i,j)%control_rod%lowered_to_level=getControlRodLoweredToLevel(sim_config, i, j)
        else if (sim_config%channel_layout_config(i,j) == CONFIG_EMPTY) then
          reactor_core(i,j)%type_of_channel=EMPTY_CHANNEL
        else if (sim_config%channel_layout_config(i,j) == CONFIG_NEUTRON_GENERATOR) then
          reactor_core(i,j)%type_of_channel=NEUTRON_GENERATOR
          ! Half a gram per cm in height
          reactor_core(i,j)%neutron_generator%weight=sim_config%size_z*100 * NEUTRON_GENERATOR_WEIGHT_PER_CM
        end if
      end do
      do j=sim_config%num_channel_configs(i), sim_config%channels_y-1
        reactor_core(i,j)%type_of_channel=EMPTY_CHANNEL
        reactor_core(i,j)%x_centre=centre_x
        reactor_core(i,j)%y_centre=((j*CHANNEL_SIZE) + (CHANNEL_SIZE/2))
      end do
    end do
  end subroutine initialiseReactorCore

  ! Initialises the neutron storage data structures at the start of the simulation so that all
  ! neutrons are inactive. For performance we hold an array of indexes, each of which represents
  ! empty slots in the neutron array which can be used. The currentNeutronIndex variable points
  ! to the currnet end of the list, and when adding a neutron the value in that location of neutron_index
  ! is read and currentNeutronIndex decremented. When deactivating a neutron the index of the newly freed
  ! location is added to the next element of neutron_index and currentNeutronIndex is incremented
  subroutine initialiseNeutrons(sim_config)
    type(simulation_configuration_type), intent(in) :: sim_config

    integer :: i

    allocate(neutrons(sim_config%max_neutrons))
    allocate(neutron_index(sim_config%max_neutrons))

    do i=1, sim_config%max_neutrons
      neutron_index(i)=i
      neutrons(i)%active=.false.
    end do
    currentNeutronIndex=sim_config%max_neutrons+1
  end subroutine initialiseNeutrons

  ! For a control rod channel will return the absolute z height position that this has been lowered to based upon
  ! the percentage insertion that was configured
  real(kind=8) function getControlRodLoweredToLevel(sim_config, channel_x, channel_y)
    type(simulation_configuration_type), intent(in) :: sim_config
    integer, intent(in) :: channel_x, channel_y

    integer rodConfigurationIndex

    rodConfigurationIndex=findControlRodConfiguration(sim_config, channel_x, channel_y)
    if (rodConfigurationIndex .lt. 0) then
      print *, "Expected control rod configuration for x=", channel_x, " y=", channel_y, &
        " but none can be found"
      call exit(-1)
    end if
    getControlRodLoweredToLevel=sim_config%size_z*(&
      sim_config%control_rod_configurations(rodConfigurationIndex)%percentage  / 100.0)
  end function getControlRodLoweredToLevel

  ! Writes out the current state of the reactor at this timestep to a file
  subroutine writeReactorState(sim_config, timestep, filename)
    type(simulation_configuration_type), intent(in) :: sim_config
    integer, intent(in) :: timestep
    character(len=*), intent(in) :: filename

    integer*8 :: num_fissions
    real(kind=8) :: mev, joules, pc(11)
    logical :: file_exists
    integer :: i, j

    num_fissions=getTotalNumberFissions(sim_config)
    mev=getMeVFromFissions(num_fissions)
    joules=getJoulesFromMeV(mev)

    inquire(file=filename, exist=file_exists)
    if (file_exists) then
      open(unit=1, file=filename, status="old", position="append", action="write")
    else
      open(unit=1, file=filename, status="unknown", action="write")
    end if

    write(1,*) "Reactor state at time ", trim(adjustl(realToExprString((NS_AS_SEC*sim_config%dt)*timestep))), &
      " secs, ",trim(adjustl(longToString(num_fissions)))," fissions releasing ", &
      trim(adjustl(realToExprString(mev)))," MeV and ", trim(adjustl(realToExprString(joules))), " Joules"
    write(1,*) "----------------------------------------------------------------------------"
    do i=0, sim_config%channels_x-1
      do j=0, sim_config%channels_y-1
        if (reactor_core(i,j)%type_of_channel==FUEL_ASSEMBLY) then
          call getFuelAssemblyChemicalContents(reactor_core(i,j)%fuel_assembly, pc)
          write(1,*) "Fuel assembly ", trim(adjustl(intToString(i))), " ", trim(adjustl(intToString(j))), &
            " ", trim(adjustl(realToString(pc(1)))), " U235 ", trim(adjustl(realToExprString(pc(2)))), " U238 ",&
            trim(adjustl(realToExprString(pc(3)))), " Pu239 ", trim(adjustl(realToExprString(pc(4)))), " U236 ", &
            trim(adjustl(realToExprString(pc(5)))), " Ba141 ", trim(adjustl(realToExprString(pc(6)))), " Kr92 ", &
            trim(adjustl(realToExprString(pc(7)))), " Xe140 ", trim(adjustl(realToExprString(pc(8)))), " Sr94 ", &
            trim(adjustl(realToExprString(pc(9)))), " Xe134 ", trim(adjustl(realToExprString(pc(10)))), " Zr103 ", &
            trim(adjustl(realToExprString(pc(11)))), " Pu240"
        end if
      end do
    end do
    write(1,*) "==========================================================================="
    close(1)
  end subroutine writeReactorState

  ! Retrieves the quantities of atoms in a fuel assembly across all the pellets for
  ! each chemical that will be written out to the file
  subroutine getFuelAssemblyChemicalContents(fuel_rod, amounts)
    type(fuel_assembly_type), intent(in) :: fuel_rod
    real(kind=8), intent(out) :: amounts(11)

    integer :: i, j

    do i=1, 11
      amounts(i)=0
    end do

    do i=1, fuel_rod%num_pellets
      do j=1, 11
        amounts(j)=amounts(j) + fuel_rod%quantities(j,i)
      end do
    end do
  end subroutine getFuelAssemblyChemicalContents

  ! Clears out the file that we are going to write to for the reactor state, this is called at simulation
  ! startup and it will overwrite any existing contents
  subroutine clearReactorStateFile(filename)
    character(len=*), intent(in) :: filename

    open(unit=1, file=filename, status="unknown", action="write")
    close(1, status='delete')
  end subroutine clearReactorStateFile

  ! Given an x and y position in the reactor core this will locate the channel
  ! that that corresponds to
  logical function locateChannelFromPosition(x, y, sim_config, channel_x, channel_y)
    real(kind=8), intent(in) :: x, y
    type(simulation_configuration_type), intent(in) :: sim_config
    integer, intent(out) :: channel_x, channel_y

    if (x .gt. sim_config%size_x .or. x .lt. 0.0 &
        .or. y .gt. sim_config%size_y .or. y .lt. 0.0) then
      locateChannelFromPosition=.false.
    else
      channel_x=x/0.2
      channel_y=y/0.2
      locateChannelFromPosition=.true.
    end if
  end function locateChannelFromPosition

  ! Based upon the properties of each fuel assembly will return the total number of fissions
  ! that have occured across all fuel assemblies in the simulation.
  integer*8 function getTotalNumberFissions(sim_config)
    type(simulation_configuration_type), intent(in) :: sim_config

    integer*8 :: total_fissions
    integer :: i, j

    total_fissions=0

    do i=0, sim_config%channels_x-1
      do j=0, sim_config%channels_y-1
        if (reactor_core(i,j)%type_of_channel == FUEL_ASSEMBLY) then
          total_fissions=total_fissions+reactor_core(i,j)%fuel_assembly%num_fissions
        end if
      end do
    end do
    getTotalNumberFissions=total_fissions
  end function getTotalNumberFissions

  ! Determines the number of currently active neutrons in the simulation
  integer*8 function getNumberActiveNeutrons(sim_config)
    type(simulation_configuration_type), intent(in) :: sim_config

    integer*8 :: activeNeutrons, i

    activeNeutrons=0

    do i=1, sim_config%max_neutrons
      if (neutrons(i)%active) activeNeutrons=activeNeutrons+1
    end do
    getNumberActiveNeutrons=activeNeutrons
  end function getNumberActiveNeutrons

  ! Converts an integer value into a string
  character(len=12) function intToString(intVal)
    integer, intent(in) :: intVal

    write(intToString,'(I10)') intVal
  end function intToString

    ! Converts a long value into a string
  character(len=12) function longToString(intVal)
    integer(kind=8), intent(in) :: intVal

    write(longToString,'(I10)') intVal
  end function longToString

  ! Converts a real value into a string
  character(len=14) function realToString(realVal)
    real(kind=8), intent(in) :: realVal

    write(realToString,'(f10.2)') realVal
  end function realToString

    ! Converts a real value into a string
  character(len=14) function realToExprString(realVal)
    real(kind=8), intent(in) :: realVal

    write(realToExprString,'(e10.4)') realVal
  end function realToExprString

  ! Returns the time difference (in seconds) between start and end time
  real(kind=8) function getTimeDifference(start_time, end_time)
    integer, intent(in) :: start_time(8), end_time(8)

    integer startms, endms

    startms=getMSFromTime(start_time)
    endms=getMSFromTime(end_time)
    getTimeDifference=(endms-startms)/1000.0
  end function getTimeDifference

  ! Returns milliseconds absolute time from provided date_time format
  integer(kind=8) function getMSFromTime(time_values)
		integer :: time_values(8)

		getMSFromTime = ( time_values(5) )*60
	 	getMSFromTime = ( getMSFromTime + time_values(6) )*60
		getMSFromTime = ( getMSFromTime + time_values(7) )*1e3
		getMSFromTime = getMSFromTime + time_values(8)
	end function getMSFromTime
end module reactor_simulator

! Program entry point, passes configuration file name to run simulation procedure
program reactor
  use reactor_simulator, only : run_simulation
  implicit none

  character(len=32) :: config_fn, output_fn

  if (command_argument_count() .ne. 2) then
    print *, "You must provide two arguments, the configuration filename and output filename to the code"
    call exit(-1)
  end if

  call get_command_argument(1, config_fn)
  call get_command_argument(2, output_fn)

  call run_simulation(config_fn, output_fn)
end program reactor
