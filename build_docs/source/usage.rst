Usage
=====

Background
----------

The general premise of the ``battery_sizing_model`` package is a series of objects
or "blocks" that are combined to create the battery simulation. The blocks that make
up a simulation are the classes :py:class:`~battery_sizing_model.blocks.Battery`, 
:py:class:`~battery_sizing_model.blocks.DegradationModel`, :py:class:`~battery_sizing_model.blocks.ECM`,
:py:class:`~battery_sizing_model.blocks.SOCIntegrator`, and 
:py:class:`~battery_sizing_model.blocks.OCVTable`.

The blocks are instantiated with whatever attributes or parameters necessary and then are combined
into the :py:class:`~battery_sizing_model.blocks.BatterySimulator` class.

Running a Simulation
--------------------

To run a simulation a power series with some constant timestep must be provided
and the simulation blocks must be created with their necessary parameters. An example
is shown below:

.. code-block:: python

    # import the package
    import battery_sizing_model.blocks as blocks

    # declare the timestep
    dt = 60

    # create an equivalent circuit block
    ecm = blocks.ECM(
        num_rc_pairs=3, 
        R0=0.00851,  
        R=[1.25430794e-03,2.75299761e-03,2.17698493e-03], 
        C=[1.65306594e+02,1.57444756e-01,1.42594864e+00], 
        timestep=dt
        )

    # create an open circuit voltage lookup table block
    ocv = blocks.OCVTable()

    # create a state of charge integrator block
    integrator = blocks.SOCIntegrator()

    # create a SOH estimation model block
    deg_model = blocks.DegradationModel()

    # create a Battery block
    battery = blocks.Battery(
        cell_chemistry="NCA", 
        pack_capacity=200e3, # Wh 
        cell_charge_energy_eff=0.9, 
        cell_discharge_energy_eff=0.9
        )

    # add the blocks to simulation instance
    sim = blocks.BatterySimulator(
        ecm=ecm,
        ocv=ocv,
        integrator=integrator,
        deg_model=deg_model,
        battery=battery,
        soc=0.5
        )

    # run the simulation assuming that there is a power
    # series named ``power``
    sim.run(power=power, timestep=dt)

A similar example is provided on GitHub in the `examples section <https://github.com/seanbuchanan-eng/battery_sizing_model/blob/main/examples/battery_sim_example.ipynb>`_
along with example data files for what the expected power series would be. For more information on attributes that can be set in
the objects see the :doc:`api`.

Gotchas
-------

There are a few assumptions that were made during the development of this simulation model that need to be discussed.

First, the :py:class:`~battery_sizing_model.blocks.Battery` class converts the power signal to power series to a cell 
level current in Amps using the nominal cell voltage of the pack. The significance of this is that the currents at high 
voltages will be higher than real and the currents at low voltages will be lower than real. It's a reasonable assumption
to make; however, it needs to be understood when analyzing results like over or under voltage alarms.

Second, the degradation model in the :py:class:`~battery_sizing_model.blocks.DegradationModel` class is not strictly 
monotonically decreasing. If the input power series has long durations of differing mean DOD or mean SOC then the model
may see a significant recovery in SOH that likely wouldn't happen in real life. It is up to the user to recognize these
situations and think critically about them. Improving this aspect of the model is something to be targeted in future work. 
 