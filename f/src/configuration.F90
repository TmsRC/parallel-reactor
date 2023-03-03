module configuration
  implicit none

  ! Maximum number of chemicals that will be configured
  integer, parameter :: NUM_CONFIG_CHEMICALS=14
  ! Maximum number of reactor core channels in a row
  integer, parameter :: MAX_CHANNEL_COLUMNS=1000
  ! Maximum number of control rods
  integer, parameter :: MAX_CONTROL_RODS=1000

  ! Type of moderator material
  integer, parameter :: WATER_MOD_TYPE_CONFIG=0, DEUTERIUM_MOD_TYPE_CONFIG=1, &
                        GRAPHITE_MOD_TYPE_CONFIG=2, NONE_MOD_TYPE_CONFIG=3

  ! Types of channel as described in configuration file
  integer, parameter :: CONFIG_EMPTY=0, CONFIG_FUEL_ASSEMBLY=1, CONFIG_CONTROL_ROD=2, &
      CONFIG_NEUTRON_GENERATOR=4, CONFIG_MODERATOR=5

  ! Different chemicals that can be in the configuration file
  integer, parameter :: U235_CONFIG=0, U238_CONFIG=1, Pu239_CONFIG=2, U236_CONFIG=3, &
                        Ba141_CONFIG=4, Kr92_CONFIG=5, Xe140_CONFIG=6, Sr94_CONFIG=7, &
                        Xe134_CONFIG=8, Zr103_CONFIG=9, Pu240_CONFIG=10, H2O_CONFIG=11, &
                        D2O_CONFIG=12, C6_CONFIG=13, UNKNOWN_CHEMICAL_CONFIG=14

  ! Holds configuration of control rod
  type control_rod_config_type
    integer :: x_channel, y_channel, percentage
  end type control_rod_config_type

  ! Overall configuration of the simulation
  type simulation_configuration_type
    integer :: num_timesteps, dt, size_x, size_y, size_z, display_progess_frequency, &
              write_reactor_state_frequency, channels_x, channels_y, moderator_weight, &
              collision_prob_multiplyer, num_ctrl_rod_configurations
    integer*8 :: max_neutrons
    integer :: moderator_type
    integer, dimension(:), allocatable:: num_channel_configs, fuel_makeup_percentage
    integer, dimension(:,:), allocatable :: channel_layout_config
    type(control_rod_config_type), dimension(:), allocatable :: control_rod_configurations
  end type simulation_configuration_type
contains

  ! A simple configuration file reader, I don't think you will need to change this (but feel free if you want to!)
  ! It will parse the configuration file and set the appropriate configuration points that will then feed into the simulation
  ! setup. It is somewhat limited in its flexibility and you need to be somewhat careful about the configuration file format,
  ! but is fine for our purposes
  subroutine parseConfiguration(config_file, sim_config)
    character(len=*), intent(in) :: config_file
    type(simulation_configuration_type), intent(out) :: sim_config

    character(len=1000) :: buffer
    integer :: file_status, entity_number, i, equals_loc
    logical :: continue_parsing, has_value

    call initialiseSimulationConfiguration(sim_config)

    open(unit=1, file=config_file, status='old', access='sequential', form='formatted', &
      action='read', iostat=file_status)

    continue_parsing=.true.
    do while (continue_parsing)
      read(1, '(A)', iostat=file_status) buffer
      if (is_iostat_end(file_status)) then
        continue_parsing=.false.
      else
        buffer=adjustl(buffer)
        ! If this starts with a hash (comment) then proceed to next line
        if (buffer(1:1) .eq. "#") cycle
        if (len_trim(buffer) .gt. 0) then
          if (hasValue(buffer)) then
            has_value=getIntValue(buffer, "NUM_TIMESTEPS", sim_config%num_timesteps)
            has_value=getIntValue(buffer, "WRITE_REACTOR_STATE_FREQUENCY", sim_config%write_reactor_state_frequency)
            has_value=getIntValue(buffer, "DISPLAY_PROGRESS_FREQUENCY", sim_config%display_progess_frequency)
            has_value=getIntValue(buffer, "DT", sim_config%dt)
            has_value=getIntValue(buffer, "SIZE_X", sim_config%size_x)
            has_value=getIntValue(buffer, "SIZE_Y", sim_config%size_y)
            has_value=getIntValue(buffer, "SIZE_Z", sim_config%size_z)
            has_value=getIntValue(buffer, "MOD_WEIGHT", sim_config%moderator_weight)
            has_value=getIntValue(buffer, "COLLISION_PROB_MULTIPLYER", sim_config%collision_prob_multiplyer)
            has_value=getLongIntValue(buffer, "MAX_NEUTRONS", sim_config%max_neutrons)

            if (hasKey(buffer, "CHANNELROW_", .true.)) then
              entity_number=getEntityNumber(buffer)
              if (entity_number .ge. 0) then
                call setChannelRowConfig(sim_config, buffer, entity_number)
              else
                print *, "Ignoring body configuration line '", trim(adjustl(buffer)),&
                      "' as this is malformed and can not extract channel row number"
              end if
            else if (hasKey(buffer, "CONTROLROD_", .true.)) then
              call setControlRodConfiguration(sim_config, buffer)
            else if (hasKey(buffer, "FUEL_", .true.)) then
              sim_config%fuel_makeup_percentage(getChemical(buffer))=getIntPercentValue(buffer)
            else if (hasKey(buffer, "MODERATOR", .true.)) then
              equals_loc=index(buffer, "=")
              if (buffer(equals_loc+1:len_trim(buffer)) == "WATER") then
                sim_config%moderator_type=WATER_MOD_TYPE_CONFIG
              else if (buffer(equals_loc+1:len_trim(buffer)) == "DEUTERIUM") then
                sim_config%moderator_type=DEUTERIUM_MOD_TYPE_CONFIG
              else if (buffer(equals_loc+1:len_trim(buffer)) == "GRAPHITE") then
                sim_config%moderator_type=GRAPHITE_MOD_TYPE_CONFIG
              else
                sim_config%moderator_type=NONE_MOD_TYPE_CONFIG
              end if
            end if
          end if
        end if
      end if
    end do
    close(1)

    sim_config%channels_x=(sim_config%size_x * 100) / 20
    sim_config%channels_y=(sim_config%size_y * 100) / 20
    call checkConfiguration(sim_config)
  end subroutine parseConfiguration

  ! Finds the control rod configuration index for a specific channel in x and y
  integer function findControlRodConfiguration(sim_config, channel_x, channel_y)
    type(simulation_configuration_type), intent(in) :: sim_config
    integer, intent(in) :: channel_x, channel_y

    integer :: i

    do i=0, sim_config%num_ctrl_rod_configurations-1
      if (sim_config%control_rod_configurations(i)%x_channel == channel_x .and. &
        sim_config%control_rod_configurations(i)%y_channel == channel_y) then
        findControlRodConfiguration=i
        return
      end if
    end do
    findControlRodConfiguration=-1
  end function findControlRodConfiguration

  ! Sets the control rod at channel x and y as per _x_y to be inserted by a specified percentage
  subroutine setControlRodConfiguration(sim_config, buffer)
    type(simulation_configuration_type), intent(inout) :: sim_config
    character(len=*), intent(in) :: buffer

    integer :: underScoreLocation, secondUnderScoreLocation, equalsLocation, &
      x_channel, y_channel, rod_index

    underScoreLocation=index(buffer, "_")
    if (underScoreLocation .gt. 0) then
      secondUnderScoreLocation=index(buffer(underScoreLocation+1:len_trim(buffer)), "_")
      if (secondUnderScoreLocation .gt. 0) then
        secondUnderScoreLocation=secondUnderScoreLocation+underScoreLocation
        equalsLocation=index(buffer(secondUnderScoreLocation+1:len_trim(buffer)), "=")
        equalsLocation=equalsLocation+secondUnderScoreLocation
        x_channel=toInt(buffer(underScoreLocation+1:secondUnderScoreLocation-1))
        y_channel=toInt(buffer(secondUnderScoreLocation+1:equalsLocation-1))
        rod_index=sim_config%num_ctrl_rod_configurations
        sim_config%control_rod_configurations(rod_index)%x_channel=x_channel
        sim_config%control_rod_configurations(rod_index)%y_channel=y_channel
        sim_config%control_rod_configurations(rod_index)%percentage=getIntPercentValue(buffer)
        sim_config%num_ctrl_rod_configurations=sim_config%num_ctrl_rod_configurations+1
      else
        print *, "Control rod percentage configuration malformed '", trim(adjustl(buffer)), "'"
        call exit(-1)
      end if
    else
      print *, "Control rod percentage configuration malformed '", trim(adjustl(buffer)), "'"
      call exit(-1)
    end if
  end subroutine setControlRodConfiguration

  ! Once configuration has been read then check it for validity before returning back to the main code, this
  ! helps ensure that errornous values, or combinations of them, are not used in simulation
  subroutine checkConfiguration(sim_config)
    type(simulation_configuration_type), intent(inout) :: sim_config

    integer :: i, x_channel, y_channel, channel_pc

    do i=0, NUM_CONFIG_CHEMICALS-1
      ! Check percentage of each chemical in fuel is a valid percentage value
      if (sim_config%fuel_makeup_percentage(i) .gt. 100 .or. &
        sim_config%fuel_makeup_percentage(i) .lt. 0) then
        print *, "Reactor fuel provided at percentage ", sim_config%fuel_makeup_percentage(i), &
          " but this is an incorrect percentage out of 100"
        call exit(-1)
      end if
    end do

    do i=0, sim_config%num_ctrl_rod_configurations-1
      ! Check control rod configuration matches a control rod channel and percentage is valid
      x_channel=sim_config%control_rod_configurations(i)%x_channel
      y_channel=sim_config%control_rod_configurations(i)%y_channel
      channel_pc=sim_config%control_rod_configurations(i)%percentage
      if (channel_pc .gt. 100 .or. channel_pc .lt. 0) then
        print *, "Control rod at x=", x_channel, " y=", y_channel, " specified at percentage ", &
          channel_pc, " but this is an incorrect percentage out of 100"
        call exit(-1)
      end if
      if (sim_config%num_channel_configs(x_channel) .gt. y_channel) then
        if (sim_config%channel_layout_config(x_channel, y_channel) .ne. CONFIG_CONTROL_ROD) then
          print *, "Control rod configuration provided for channel x=", x_channel," y=", y_channel, &
            " but this type is not a control rod"
          call exit(-1)
        end if
      else
        print *, "Control rod configuration provided for channel x=", x_channel, " y=", y_channel, &
          " but this is an empty channel"
        call exit(-1)
      end if
    end do

    do i=0, MAX_CHANNEL_COLUMNS-1
      ! Check channels that are configured are sensible
      if (sim_config%num_channel_configs(i) .gt. sim_config%channels_y) then
        print *, "There are ", sim_config%channels_y," reactor channel columns, but row ", i, &
          " has ", sim_config%num_channel_configs(i)
        call exit(-1)
      end if
    end do
  end subroutine checkConfiguration

  ! Initialises the configuration data structure for this simulation, adds in some defaults
  ! if they are not specified and marks each body as inactive
  subroutine initialiseSimulationConfiguration(sim_config)
    type(simulation_configuration_type), intent(out) :: sim_config

    integer :: i

    sim_config%dt=10
    sim_config%num_timesteps=1000
    sim_config%num_ctrl_rod_configurations=0
    sim_config%write_reactor_state_frequency=0
    sim_config%display_progess_frequency=0
    ! Note that we allocate all these arrays as zero indexed because that is the
    ! starting index used in the configuration file, and therefore simpler to keep
    ! the same as the C code here for this.
    allocate(sim_config%channel_layout_config(0:MAX_CHANNEL_COLUMNS-1, 0:MAX_CHANNEL_COLUMNS-1))
    allocate(sim_config%control_rod_configurations(0:MAX_CONTROL_RODS-1))
    allocate(sim_config%num_channel_configs(0:MAX_CHANNEL_COLUMNS-1))
    allocate(sim_config%fuel_makeup_percentage(0:NUM_CONFIG_CHEMICALS-1))
    do i=0, MAX_CHANNEL_COLUMNS-1
      sim_config%num_channel_configs(i)=0
    end do
    do i=0, NUM_CONFIG_CHEMICALS-1
      sim_config%fuel_makeup_percentage(i)=0
    end do
  end subroutine initialiseSimulationConfiguration

  ! Retrieves the chemical index that has been added after the underscore,
  !  for instance FUEL_U235 will return 0
  integer function getChemical(buffer)
    character(len=*), intent(in) :: buffer

    integer :: underScoreLocation, equalsLocation, from_idx, to_idx

    underScoreLocation=index(buffer, "_")
    if (underScoreLocation .gt. 0) then
      equalsLocation=index(buffer(underScoreLocation+1:len_trim(buffer)), "=")
      if (equalsLocation .gt. 0) then
        from_idx=underScoreLocation+1
        to_idx=underScoreLocation+equalsLocation-1
        if (buffer(from_idx:to_idx) == "U235") then
          getChemical=U235_CONFIG
        else if (buffer(from_idx:to_idx) == "U238") then
          getChemical=U238_CONFIG
        else if (buffer(from_idx:to_idx) == "Pu239") then
          getChemical=Pu239_CONFIG
        else if (buffer(from_idx:to_idx) == "U236") then
          getChemical=U236_CONFIG
        else if (buffer(from_idx:to_idx) == "Ba141") then
          getChemical=Ba141_CONFIG
        else if (buffer(from_idx:to_idx) == "Kr92") then
          getChemical=Kr92_CONFIG
        else if (buffer(from_idx:to_idx) == "Xe140") then
          getChemical=Xe140_CONFIG
        else if (buffer(from_idx:to_idx) == "Sr94") then
          getChemical=Sr94_CONFIG
        else if (buffer(from_idx:to_idx) == "Xe134") then
          getChemical=Xe134_CONFIG
        else if (buffer(from_idx:to_idx) == "Zr103") then
          getChemical=Zr103_CONFIG
        else if (buffer(from_idx:to_idx) == "Pu240") then
          getChemical=Pu240_CONFIG
        else
          getChemical=UNKNOWN_CHEMICAL_CONFIG
        end if
      else
        getChemical=UNKNOWN_CHEMICAL_CONFIG
      end if
    else
      getChemical=UNKNOWN_CHEMICAL_CONFIG
    end if
  end function getChemical

  ! Based upon a source string will figure out the reactor core's channels configurations
  ! and return this, along with the number of options specified. For instance
  ! "FUEL,CONTROL,MODERATOR,CONTROL,FUEL" would return enumerations matching those
  ! strings and five as the number of options
  subroutine setChannelRowConfig(sim_config, buffer, row)
    type(simulation_configuration_type), intent(inout) :: sim_config
    character(len=*), intent(in) :: buffer
    integer, intent(in) :: row

    integer :: number_options_specified, i, comma_location, start_loc, &
              equals_loc
    character(len=100) :: option_str

    equals_loc=index(buffer, "=")
    if (equals_loc .gt. 0) then
      start_loc=equals_loc+1
      number_options_specified=countNumOccurances(buffer(start_loc:len_trim(buffer)), ',')
      sim_config%num_channel_configs(row)=number_options_specified+1
      do i=0, number_options_specified
        if (i .lt. number_options_specified) then
          comma_location=index(buffer(start_loc:len_trim(buffer)), ",")-1
        else
          comma_location=len_trim(buffer)
        end if
        option_str=buffer(start_loc:start_loc+comma_location-1)
        if (option_str == "FUEL") then
          sim_config%channel_layout_config(row, i)=CONFIG_FUEL_ASSEMBLY
        else if (option_str == "CONTROL") then
          sim_config%channel_layout_config(row, i)=CONFIG_CONTROL_ROD
        else if (option_str == "NEUTRON") then
          sim_config%channel_layout_config(row, i)=CONFIG_NEUTRON_GENERATOR
        else if (option_str == "MODERATOR") then
          sim_config%channel_layout_config(row, i)=CONFIG_MODERATOR
        else if (option_str == "EMPTY") then
          sim_config%channel_layout_config(row, i)=CONFIG_EMPTY
        end if
        start_loc=start_loc+comma_location+1
      end do
    end if
  end subroutine setChannelRowConfig

  ! Helper to count the number of occurances of the character 'c' in
  ! a string that has been provided
  integer function countNumOccurances(buffer, c)
    character(len=*), intent(in) :: buffer
    character, intent(in) :: c

    integer :: i, occurances

    occurances=0

    do i=1, len_trim(buffer)
      if (buffer(i:i) == c) occurances=occurances+1
    end do
    countNumOccurances=occurances
  end function countNumOccurances

  ! Given a string with a key-value pair (e.g. key = value) this will extract the integer value after the equals
  ! and return this via the dummy argument. It returns true if such extraction was possible and
  ! false if not
  logical function getIntValue(buffer, configuration_key, configuration_value)
    character(len=*), intent(in) :: buffer, configuration_key
    integer, intent(out) :: configuration_value

    integer :: key_loc, equals_loc

    equals_loc=index(buffer, "=")
    if (equals_loc .gt. 0) then
      key_loc=index(buffer, configuration_key)
      if (key_loc .gt. 0) then
        configuration_value=toInt(buffer(equals_loc+1:len_trim(buffer)))
        getIntValue=.true.
        return
      end if
    end if
    getIntValue=.false.
  end function getIntValue

  ! Given a string with a key-value pair (e.g. key = value) this will extract the long integer
  ! value after the equals and return this via the dummy argument. It returns true if such extraction
  ! was possible and false if not
  logical function getLongIntValue(buffer, configuration_key, configuration_value)
    character(len=*), intent(in) :: buffer, configuration_key
    integer*8, intent(out) :: configuration_value

    integer :: key_loc, equals_loc

    equals_loc=index(buffer, "=")
    if (equals_loc .gt. 0) then
      key_loc=index(buffer, configuration_key)
      if (key_loc .gt. 0) then
        configuration_value=toLongInt(buffer(equals_loc+1:len_trim(buffer)))
        getLongIntValue=.true.
        return
      end if
    end if
    getLongIntValue=.false.
  end function getLongIntValue

  ! From a string will extract the percentage integer value after the
  ! equals and return this, or -1 if none is found
  integer function getIntPercentValue(buffer)
    character(len=*), intent(in) :: buffer

    integer :: equals_loc, pcLocation

    equals_loc=index(buffer, "=")
    if (equals_loc .gt. 0) then
      pcLocation=index(buffer, "%")
      if (pcLocation .gt. 0) then
        getIntPercentValue=toInt(buffer(equals_loc+1:pcLocation-1))
        return
      else
        getIntPercentValue=toInt(buffer(equals_loc+1:len_trim(buffer)))
        return
      end if
    end if
    getIntPercentValue=-1
  end function getIntPercentValue

  ! A helper function to parse a string with an underscore in it, this will extract the number after
  ! the underscore as we use this in the configuration file for setting numbers of bodies in the configuration
  integer function getEntityNumber(buffer)
    character(len=*), intent(in) :: buffer

    integer :: uscore_loc, equals_loc

    uscore_loc=index(buffer, "_")
    if (uscore_loc .gt. 0) then
      equals_loc=index(buffer(uscore_loc+1:), "=")
      if (equals_loc .gt. 0) then
        getEntityNumber=toInt(buffer(uscore_loc+1:uscore_loc+equals_loc-1))
        return
      end if
    end if
    print *, "Line '", trim(adjustl(buffer)), "' is malformed, can not extract entity number"
    call exit(-1)
    getEntityNumber=-1
  end function getEntityNumber

  ! Given a string with a key-value pair (e.g. key = value) this will extract the real value after the equals
  ! and return this via the dummy argument. It returns true if such extraction was possible and
  ! false if not
  logical function getRealValue(buffer, configuration_key, configuration_value)
    character(len=*), intent(in) :: buffer, configuration_key
    real(kind=8), intent(out) :: configuration_value

    integer :: key_loc, equals_loc

    equals_loc=index(buffer, "=")
    if (equals_loc .gt. 0) then
      key_loc=index(buffer, configuration_key)
      if (key_loc .gt. 0) then
        configuration_value=toReal(buffer(equals_loc+1:len_trim(buffer)))
        getRealValue=.true.
        return
      end if
    end if
    getRealValue=.false.
  end function getRealValue

  ! Determines whether a configuration value holds a key value pair (and specifically whether it has a value)
  logical function hasValue(buffer)
    character(len=*), intent(in) :: buffer

    hasValue=index(buffer, "=") .gt. 0
  end function hasValue

  ! Determines whether the current line represented by buffer has the provided
  ! configuration key, with the function returning either true or false.
  ! There is an optional restriction that can be imposed, namely whether that
  ! key starts at the beginning of the line or not (if omitted the default is false for this flag)
  logical function hasKey(buffer, configuration_key, from_start)
    character(len=*), intent(in) :: buffer, configuration_key
    logical, optional, intent(in) :: from_start

    integer :: key_loc, equals_loc
    equals_loc=index(buffer, "=")
    if (equals_loc .gt. 0) then
      key_loc=index(buffer, configuration_key)
      if (key_loc .gt. 0) then
        if (present(from_start) .and. from_start) then
          hasKey=key_loc .eq. 1
        else
          hasKey=.true.
        end if
        return
      end if
    end if
    hasKey=.false.
  end function hasKey

  ! Converts a string to an integer
  integer function toInt(str)
    character(len=*), intent(in) :: str

    read(str,*)  toInt
  end function toInt

  ! Converts a string to an integer
  integer*8 function toLongInt(str)
    character(len=*), intent(in) :: str

    read(str,*)  toLongInt
  end function toLongInt

  real(kind=8) function toReal(str)
    character(len=*), intent(in) :: str

    read(str,*)  toReal
  end function toReal
end module configuration
