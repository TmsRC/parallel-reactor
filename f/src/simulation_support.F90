module simulation_support
  implicit none

  ! Neutrons per gram of Californium 252 released per second
  integer*8, parameter :: CF252_NEUTRONS_PER_GRAM_PER_SEC=23e12
  ! Number of nanoseconds in a second
  real(kind=8), parameter :: NS_AS_SEC=1e-9
  ! Speed of light in m/s
  integer, parameter :: C=3e8
  !Constant of proportionality
  real(kind=8), parameter :: R0=1.2e-15
  ! PI
  real(kind=8), parameter :: PI=3.14159265359
  ! Planck's constant multiplied by speed of light
  real(kind=8), parameter :: HC=1239.84
  ! Size of a barn in meters squared
  real(kind=8), parameter :: BARN=1e-28
  ! The Avogadro constant
  real(kind=8), parameter :: AVOGADRO=6.022e23
  ! Amount of energy in MeV released per fission
  integer, parameter :: MeV_PER_FISSION=200
  ! Number of Joules of energy corresponding to 1 MeV
  real(kind=8), parameter :: JOULES_PER_MEV=1.6e-13
  ! The scattering cross section of hydrogen
  integer, parameter :: SCATTERING_CROSS_SECTION_HYDROGEN=20
  ! The scattering cross section of oxygen
  integer, parameter :: SCATTERING_CROSS_SECTION_OXYGEN=3.8
  ! The scattering cross section of carbon
  integer, parameter :: SCATTERING_CROSS_SECTION_CARBON=4.9
  ! Density of water in grams per cm^3
  real(kind=8), parameter :: DENSITY_WATER=1.0
  ! Density of deuterium in grams per cm^3
  real(kind=8), parameter :: DENSITY_DEUTERIUM=1.8
  ! Density of graphite in grams per cm^3
  real(kind=8), parameter :: DENSITY_GRAPHITE=2.25
  ! Probability that water will absorb neutron
  real(kind=8), parameter :: H2O_ABSORBTION_PROB=0.25
  ! Probability that deuterium will absorb neutron
  real(kind=8), parameter :: D2O_ABSORPTION_PROB=0.05
  ! Probability that graphite will absorb neutron
  real(kind=8), parameter :: C_ABSORPTION_PROB=0.01
  ! Mass in kg of one unit of atomic rest mass
  real(kind=8), parameter :: MASS_ONE_UNIT=1.6726*1e-27
  ! Newton meters per electron volt
  real(kind=8), parameter :: NEWTON_METER_PER_eV=1.60217733*1e-19

  ! Total number of chemicals that the simulation can represent
  integer, parameter :: NUM_CHEMICALS=14

  ! Type of reactor channel
  integer, parameter :: MODERATOR=0, FUEL_ASSEMBLY=1, CONTROL_ROD=2, EMPTY_CHANNEL=3, NEUTRON_GENERATOR=4
  ! Type of moderator material
  integer, parameter :: WATER=0, DEUTERIUM=1, GRAPHITE=2, NO_MODERATOR=3
  ! Chemicals that can be represented in the simulation
  integer, parameter :: U235=1, U238=2, Pu239=3, U236=4, Ba141=5, Kr92=6, Xe140=7, &
                        Sr94=8, Xe134=9, Zr103=10, Pu240=11, H2O=12, D2O=13, C6=14, UNKNOWN_CHEMICAL=15

  ! Holds the state of the moderator channel
  type moderator_type
    integer :: moderator_material
  end type moderator_type

  ! Holds the state of a fuel assembly in a channel, fuel assemblies
  ! are made up of fuel pellets (40mm by 40mm by 2mm) stacked on top of each other
  type fuel_assembly_type
    real(kind=8), dimension(:,:), allocatable :: quantities
    integer :: num_pellets
    integer*8:: num_fissions
  end type fuel_assembly_type

  ! State of a control rod in a channel
  type control_rod_type
    real(kind=8) :: lowered_to_level
  end type control_rod_type

  ! State of a neutron generator in a channel
  type neutron_generator_type
    real(kind=8) :: weight
  end type neutron_generator_type

  ! Represents an individual neutron
  type neutron_type
    integer:: x,y,z
    real(kind=8) :: pos_x, pos_y, pos_z
    real(kind=8) :: energy
    logical active
  end type neutron_type

  ! A single channel in the reactor
  type channel_type
    integer :: type_of_channel
    type(moderator_type) :: moderator
    type(fuel_assembly_type) :: fuel_assembly
    type(control_rod_type) :: control_rod
    type(neutron_generator_type) :: neutron_generator
    real(kind=8) :: x_centre, y_centre
  end type channel_type

contains

  ! Initialises the simulation support by seeding the random number generator.
  ! Random seeding is based on the example at https://gcc.gnu.org/onlinedocs/gcc-4.6.4/gfortran/RANDOM_005fSEED.html
  subroutine seedRandomNumberGenerator()
    integer :: clock, n, i
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))
    call system_clock(count=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)
    deallocate(seed)
  end subroutine seedRandomNumberGenerator

  ! Given the weight of a sample of Californium 252, will calculate how many
  ! neutrons this has produced in a number of nanoseconds
  integer*8 function getNumberNeutronsFromGenerator(weight, dt)
    real(kind=8), intent(in) :: weight
    integer, intent(in) :: dt

    getNumberNeutronsFromGenerator=CF252_NEUTRONS_PER_GRAM_PER_SEC * &
      weight * (dt * NS_AS_SEC)
  end function getNumberNeutronsFromGenerator

  ! Converts energy in MeV with a specific atomic rest mass to velocity
  ! in meters per second. A neutron always has a rest mass of one
  real(kind=8) function MeVToVelocity(MeV, rest_mass)
    real(kind=8), intent(in) :: MeV
    integer, intent(in) :: rest_mass

    real(kind=8) :: mass, ke, vvc

    mass=rest_mass*MASS_ONE_UNIT
    ke=MeV*1e6*NEWTON_METER_PER_eV
    vvc=sqrt(1-1/(ke/(mass*C*C)-(-1))**2)
    MeVToVelocity=vvc * C
  end function MeVToVelocity

  ! Will determine and handle if a neutron has collided with the moderator, returning true
  ! if the neutron is no longer involved in the simulation (i.e. it has been absorbed) or
  ! false if it is still present. On scattering the energy of the neutron is updated.
  ! Note that for the absorption we use predefined probabilities, that are fairly close to those
  ! measured experimentally rather than computing the absorption cross section directly.
  logical function determineAndHandleIfNeutronModeratorCollision(neutron, moderator_weight, moderator_type, size_z)
    type(neutron_type), intent(inout) :: neutron
    integer, intent(in) :: moderator_weight, moderator_type, size_z

    real(kind=8) :: absorption_prob, mean_log_energy_reduction
    integer :: atomic_mass

    if (determineNeutronModeratorScattering(moderator_weight, moderator_type, size_z)) then
      if (moderator_type == WATER) then
        absorption_prob=H2O_ABSORBTION_PROB
        atomic_mass=1
      end if

      if (moderator_type == DEUTERIUM) then
        absorption_prob=D2O_ABSORPTION_PROB
        atomic_mass=2
      end if

      if (moderator_type == GRAPHITE) then
        absorption_prob=C_ABSORPTION_PROB
        atomic_mass=12
      end if

      if (getRandomReal() >= absorption_prob) then
        ! Use 1 here for the mass as its the hydrogen in the H20 that does the slowing
        mean_log_energy_reduction=2/(1+(atomic_mass/3))
        neutron%energy=neutron%energy/exp(mean_log_energy_reduction)
        determineAndHandleIfNeutronModeratorCollision=.false.
        return
      else
        ! Neutron absorbed by water
        determineAndHandleIfNeutronModeratorCollision=.true.
        return
      end if
    end if
    determineAndHandleIfNeutronModeratorCollision=.false.
  end function determineAndHandleIfNeutronModeratorCollision

  ! Determines and handles if a neutron has collided with an atom of fuel, returning true if
  ! so (the neutron will be deactivated from the simulation) and false otherwise. This handles both
  ! the U235 and Pu239 fuels (it is assumed that no other fuels fission).
  logical function determineAndHandleIfNeutronFuelCollision(MeV, reactor_channel, fuel_pellet, collision_prob_multiplyer)
    real(kind=8), intent(in) :: MeV
    type(channel_type), intent(inout) :: reactor_channel
    integer, intent(in) :: fuel_pellet, collision_prob_multiplyer

    real(kind=8) :: deBroglieWavelength

    deBroglieWavelength=calculateDeBroglieWavelength(MeV, 1)
    ! Calculate the cross section in Barns for each chemical and determine if a collision occured
    ! if so add the neutron to the atom it collided with.
    if (reactor_channel%fuel_assembly%quantities(U235, fuel_pellet) .gt. 0) then
      if (determineNeutronAbsorbtionByFuel(deBroglieWavelength, reactor_channel, U235, &
            fuel_pellet, collision_prob_multiplyer)) then
        reactor_channel%fuel_assembly%quantities(U235, fuel_pellet) = &
          reactor_channel%fuel_assembly%quantities(U235, fuel_pellet)-1
        reactor_channel%fuel_assembly%quantities(U236, fuel_pellet) = &
          reactor_channel%fuel_assembly%quantities(U236, fuel_pellet)+1
        determineAndHandleIfNeutronFuelCollision=.true.
        return
      end if
    end if
    if (reactor_channel%fuel_assembly%quantities(Pu239, fuel_pellet) .gt. 0) then
      if (determineNeutronAbsorbtionByFuel(deBroglieWavelength, reactor_channel, Pu239, &
            fuel_pellet, collision_prob_multiplyer)) then
        reactor_channel%fuel_assembly%quantities(Pu239, fuel_pellet) = &
          reactor_channel%fuel_assembly%quantities(Pu239, fuel_pellet)-1
        reactor_channel%fuel_assembly%quantities(Pu240, fuel_pellet) = &
          reactor_channel%fuel_assembly%quantities(Pu240, fuel_pellet)+1
        determineAndHandleIfNeutronFuelCollision=.true.
        return
      end if
    end if
    determineAndHandleIfNeutronFuelCollision=.false.
  end function determineAndHandleIfNeutronFuelCollision

  ! Handles the fission of Pu240, which will either split into Xenon and
  ! Zirconium, releasing 3 neutrons, or go back to Pu239 releasing
  ! one neutron.
  integer function fissionPu240(reactor_channel, fuel_pellet)
    type(channel_type), intent(inout) :: reactor_channel
    integer, intent(in) :: fuel_pellet

    if (getRandomReal() .lt. 0.73) then
      reactor_channel%fuel_assembly%quantities(Xe134, fuel_pellet)= &
        reactor_channel%fuel_assembly%quantities(Xe134, fuel_pellet)+1
      reactor_channel%fuel_assembly%quantities(Zr103, fuel_pellet)= &
        reactor_channel%fuel_assembly%quantities(Zr103, fuel_pellet)+1
      reactor_channel%fuel_assembly%quantities(Pu240, fuel_pellet)= &
        reactor_channel%fuel_assembly%quantities(Pu240, fuel_pellet)-1
      fissionPu240=3
    else
    reactor_channel%fuel_assembly%quantities(Pu239, fuel_pellet)= &
        reactor_channel%fuel_assembly%quantities(Pu239, fuel_pellet)+1
      reactor_channel%fuel_assembly%quantities(Pu240, fuel_pellet)= &
        reactor_channel%fuel_assembly%quantities(Pu240, fuel_pellet)-1
      fissionPu240=1
    end if
  end function fissionPu240

  ! Handles the fission of U236, which will either split into BariumgetNumberAtomsInFuelPellet
  ! and Krypton releasing 3 neutrons, or Xenon and Strontium releasing
  ! 2 neutrons.
  integer function fissionU236(reactor_channel, fuel_pellet)
    type(channel_type), intent(inout) :: reactor_channel
    integer, intent(in) :: fuel_pellet

    if (getRandomReal() .lt. 0.85) then
      reactor_channel%fuel_assembly%quantities(Ba141, fuel_pellet)= &
        reactor_channel%fuel_assembly%quantities(Ba141, fuel_pellet)+1
      reactor_channel%fuel_assembly%quantities(Kr92, fuel_pellet)= &
        reactor_channel%fuel_assembly%quantities(Kr92, fuel_pellet)+1
      reactor_channel%fuel_assembly%quantities(U236, fuel_pellet)= &
        reactor_channel%fuel_assembly%quantities(U236, fuel_pellet)-1
      fissionU236=3
    else
      reactor_channel%fuel_assembly%quantities(Xe140, fuel_pellet)= &
        reactor_channel%fuel_assembly%quantities(Xe140, fuel_pellet)+1
      reactor_channel%fuel_assembly%quantities(Sr94, fuel_pellet)= &
        reactor_channel%fuel_assembly%quantities(Sr94, fuel_pellet)+1
      reactor_channel%fuel_assembly%quantities(U236, fuel_pellet)= &
        reactor_channel%fuel_assembly%quantities(U236, fuel_pellet)-1
      fissionU236=2
    end if
  end function fissionU236

  ! Initialises a neutron based on it's location in the reactor (the channel
  ! it is in, i.e. the neutron generator's channel, and it's height via z).
  subroutine initialiseNeutron(neutron, reactor_channel, z)
    type(neutron_type), intent(inout) :: neutron
    type(channel_type), intent(inout) :: reactor_channel
    real(kind=8), intent(in) :: z

    neutron%active=.true.
    call generateNeutronVelocityAndEnergy(neutron)
    neutron%pos_x=reactor_channel%x_centre
    neutron%pos_y=reactor_channel%y_centre
    neutron%pos_z=z
  end subroutine initialiseNeutron

  ! Will calculate the amount of energy in MeV that a number of fissions
  ! has produced
  real(kind=8) function getMeVFromFissions(num_fissions)
    integer*8, intent(in) :: num_fissions

    getMeVFromFissions=num_fissions*MeV_PER_FISSION
  end function getMeVFromFissions

  ! Converts Mev to Joules
  real(kind=8) function getJoulesFromMeV(MeV)
    real(kind=8), intent(in) :: MeV

    getJoulesFromMeV=JOULES_PER_MEV*MeV
  end function getJoulesFromMeV

  ! Gets the number of atoms in a gram of a specific chemical
  real(kind=8) function getAtomsPerGram(chemical)
    integer, intent(in) :: chemical

    integer*8 :: restMass

    restMass=getAtomicRestMass(chemical)
    getAtomsPerGram=AVOGADRO/restMass
  end function getAtomsPerGram

  ! Determines whether a neutron scatters from the atom of the moderator based on its chemical composition,
  ! the size of the moderator and its weight.
  !
  ! This is a simplification, as the Barn unit is metres squared
  ! whereas we assume all the atoms of the moderator are in that 2D plane, however realistically the neutron
  ! is moving through the moderator so will likely come into contact with them in the other dimension anyway
  ! in a timestep, so that is likely reasonable and it keeps things simpler too
  logical function determineNeutronModeratorScattering(moderator_weight, moderator_type, size_z)
    integer, intent(in) :: moderator_weight, moderator_type, size_z

    integer :: moderator_chemical
    real(kind=8) :: cs, moderator_barns, size_reduction, num_atoms_reduced, num_atoms

    if (moderator_type == WATER) moderator_chemical=H2O
    if (moderator_type == DEUTERIUM) moderator_chemical=D2O
    if (moderator_type == GRAPHITE) moderator_chemical=C6
    cs=getScatteringCrossSection(moderator_chemical)

    ! Number of barns of moderator height
    moderator_barns=size_z/BARN
    ! How many times larger the moderator is than the scattering cs
    size_reduction=moderator_barns/cs

    ! Get the number of atoms in the moderator
    num_atoms=getAtomsPerGram(moderator_chemical) * moderator_weight
    ! Scale the number of atoms down by the size reduction to then see if there is likely one in that square
    num_atoms_reduced=num_atoms/size_reduction
    determineNeutronModeratorScattering=num_atoms_reduced .ge. getRandomReal()
  end function determineNeutronModeratorScattering

  ! Gets the scattering cross section of a specific chemical, here we are using
  ! cross section constants that assume 1eV of energy of the neutron, which is
  ! a simplification but it keeps the maths simpler
  real(kind=8) function getScatteringCrossSection(chemical)
    integer, intent(in) :: chemical

    if (chemical == H2O) then
      getScatteringCrossSection=getScatteringCrossSectionOfWater()
    else if (chemical == D2O) then
      getScatteringCrossSection=getScatteringCrossSectionOfDeuterium()
    else if (chemical == C6) then
      getScatteringCrossSection=getScatteringCrossSectionOfGraphite()
    end if
  end function getScatteringCrossSection

  ! Get the scattering cross section of water (H2O)
  real(kind=8) function getScatteringCrossSectionOfWater()
    integer :: rest_mass
    real(kind=8) :: n

    rest_mass=getAtomicRestMass(H2O)
    n=(AVOGADRO/rest_mass) * DENSITY_WATER
    getScatteringCrossSectionOfWater=(2*n*SCATTERING_CROSS_SECTION_HYDROGEN) &
      + (n*SCATTERING_CROSS_SECTION_OXYGEN)
  end function getScatteringCrossSectionOfWater

  ! Get the scattering cross section of duterium (D2O)
  real(kind=8) function getScatteringCrossSectionOfDeuterium()
    integer :: rest_mass
    real(kind=8) :: n

    rest_mass=getAtomicRestMass(D2O)
    n=(AVOGADRO/rest_mass) * DENSITY_DEUTERIUM
    getScatteringCrossSectionOfDeuterium=(2*n*SCATTERING_CROSS_SECTION_HYDROGEN) &
      + (n*SCATTERING_CROSS_SECTION_OXYGEN)
  end function getScatteringCrossSectionOfDeuterium

  ! Get the scattering cross section of graphite (C6)
  real(kind=8) function getScatteringCrossSectionOfGraphite()
    integer :: rest_mass
    real(kind=8) :: n

    rest_mass=getAtomicRestMass(C6)
    n=(AVOGADRO/rest_mass) * DENSITY_GRAPHITE
    getScatteringCrossSectionOfGraphite=n*SCATTERING_CROSS_SECTION_CARBON
  end function getScatteringCrossSectionOfGraphite

  ! For a given neutron this will generate a random amount of energy in
  ! MeV and the velocity components (which add up to 100 across the
  ! three dimensions) that the velocity will be shared over.
  subroutine generateNeutronVelocityAndEnergy(neutron)
    type(neutron_type), intent(inout) :: neutron

    integer :: running_sum, x, y, z

    running_sum=0
    do while (running_sum == 0)
      x=getRandomInt(100)
      y=getRandomInt(100)
      z=getRandomInt(100)
      running_sum = x+y+z
    end do
    neutron%x=x*(100/running_sum)
    if (getRandomInt(2) == 1) neutron%x=-neutron%x
    neutron%y=y*(100/running_sum)
    if (getRandomInt(2) == 1) neutron%y=-neutron%y
    neutron%z=z*(100/running_sum)
    if (getRandomInt(2) == 1) neutron%z=-neutron%z
    ! Random energy between 0 and 10 MeV
    neutron%energy=getRandomReal() * 10
  end subroutine generateNeutronVelocityAndEnergy

  ! Determines whether a neutron has been absorbed by an atom of nuclear fuel within
  ! a specific fuel pellet of a fuel rod.
  !
  ! This is a simplification, as the Barn unit is metres squared
  ! whereas we assume all the atoms of the fuel pellet are in that 2D plane, however
  ! realistically the neutron is moving through the fuel channel so will likely come
  ! into contact with them in the other dimension anyway
  ! in a timestep, so that is likely reasonable and it keeps things simpler too
  logical function determineNeutronAbsorbtionByFuel(deBroglieWavelength, reactor_channel, &
      chemical, fuel_pellet, collision_prob_multiplyer)
    real(kind=8), intent(in) :: deBroglieWavelength
    type(channel_type), intent(inout) :: reactor_channel
    integer, intent(in) :: chemical, fuel_pellet, collision_prob_multiplyer

    real(kind=8) :: cs, pellet_b, size_reduction, num_atoms_reduced, num_atoms

    ! First get the cross section in Barns for the chemical
    cs=getCrossSectionInBarns(deBroglieWavelength, chemical)

    ! Determine the size of the fuel pellet in barns
    pellet_b=0.0016/BARN
    ! How much smaller the cross section is than the pellet
    size_reduction=pellet_b/cs

    ! Now grab out the number of fuel atoms in the fuel pellet and scale this down by the cross section ratio
    num_atoms=getNumberAtomsInFuelPellet(reactor_channel, chemical, fuel_pellet)
    num_atoms_reduced=num_atoms/size_reduction
    num_atoms_reduced = num_atoms_reduced * collision_prob_multiplyer
    ! Now determine randomly if there is an atom in the cross section which will result in absorption
    ! We scale the number of atoms reduced here by 1000000 to improve the probability as the Fortran
    ! random number generator doesn't tend to go down that low
    if (num_atoms_reduced .ge. getRandomReal()) then
      determineNeutronAbsorbtionByFuel=.true.
    else
      determineNeutronAbsorbtionByFuel=.false.
    end if
  end function determineNeutronAbsorbtionByFuel

  ! Retrieves the number of atoms of a chemical in a specific pellet of the fuel
  ! held in the given reactor channel
  real(kind=8) function getNumberAtomsInFuelPellet(reactor_channel, chemical, fuel_pellet)
    type(channel_type), intent(inout) :: reactor_channel
    integer, intent(in) :: chemical, fuel_pellet

    getNumberAtomsInFuelPellet=reactor_channel%fuel_assembly%quantities(chemical, fuel_pellet)
  end function getNumberAtomsInFuelPellet

  ! Based upon a collision with a chemical, will calculate the cross section in Barnes from the
  ! deBroglie wavelength that has been provided. 1 barn is 1e-28 meters squared.
  real(kind=8) function getCrossSectionInBarns(deBroglieWavelength, chemical)
    real(kind=8), intent(in) :: deBroglieWavelength
    integer, intent(in) :: chemical

    integer :: rest_mass
    real(kind=8) :: R, crossSection

    rest_mass=getAtomicRestMass(chemical)
    if (rest_mass == 0) then
      print *, "Rest mass of chemical determined as zero during cross section calculation"
      call exit(-1)
    end if
    R=R0*(rest_mass**0.03)
    crossSection=PI*((R+(deBroglieWavelength/(2*PI)))**2)
    getCrossSectionInBarns=crossSection/BARN
  end function getCrossSectionInBarns

  ! Returns the atomic rest mass of a chemical
  integer function getAtomicRestMass(chemical)
    integer, intent(in) :: chemical

    if (chemical == U235) then
      getAtomicRestMass=235
    else if (chemical == U238) then
      getAtomicRestMass=238
    else if (chemical == Pu239) then
      getAtomicRestMass=239
    else if (chemical == U236) then
      getAtomicRestMass=236
    else if (chemical == Ba141) then
      getAtomicRestMass=141
    else if (chemical == Kr92) then
      getAtomicRestMass=92
    else if (chemical == Xe140) then
      getAtomicRestMass=140
    else if (chemical == Sr94) then
      getAtomicRestMass=94
    else if (chemical == Xe134) then
      getAtomicRestMass=134
    else if (chemical == Zr103) then
      getAtomicRestMass=103
    else if (chemical == Pu240) then
      getAtomicRestMass=240
    else if (chemical == H2O) then
      getAtomicRestMass=18
    else if (chemical == D2O) then
      getAtomicRestMass=19
    else if (chemical == C6) then
      getAtomicRestMass=12
    else
      getAtomicRestMass=0
    end if
  end function getAtomicRestMass

  ! From the energy (in MeV) and atomic rest mass will calculate the
  ! deBroglie wavelength
  real(kind=8) function calculateDeBroglieWavelength(MeV, rest_mass)
    real(kind=8), intent(in) :: MeV
    integer, intent(in) :: rest_mass

    real(kind=8) :: mass, ke, ppc, pc, ddw

    mass=rest_mass*MASS_ONE_UNIT
    ke=MeV*1e6*NEWTON_METER_PER_eV
    ppc=sqrt(ke*ke-(-1)*2*ke*mass*C*C)
    pc=ppc/NEWTON_METER_PER_eV
    ddw=HC/pc
    calculateDeBroglieWavelength=ddw*NS_AS_SEC
  end function calculateDeBroglieWavelength

  ! Retrieves a random integer between 0 and max
  integer function getRandomInt(max_val)
    integer, intent(in) :: max_val

    real :: rval

    call random_number(rval)
    getRandomInt=floor(rval * max_val)
  end function getRandomInt

  ! Retrieves a random real between 0 and 1
  real(kind=8) function getRandomReal()
    real :: rval

    call random_number(rval)
    getRandomReal=real(rval, 8)
  end function getRandomReal
end module simulation_support
