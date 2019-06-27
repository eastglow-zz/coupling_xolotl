[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 120
  ny = 120
  xmin = 0
  xmax = 1200
  ymin = 0
  ymax = 1200
  #uniform_refine = 3
[]

[Variables]
  [./a]
  [../]
[]

[Executioner]
  type = Transient
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.001
    growth_factor = 1.2
    cutback_factor = 0.8
  [../]
  start_time = 0
  end_time = 7e7
  # end_time = 20000.0
[]

[Problem]
  kernel_coverage_check = false
  solve = false
[]

[MultiApps]
  [./XolotlWrapper]
    type = TransientMultiApp
    app_type = coupling_xolotlApp
    execute_on = TIMESTEP_BEGIN
    positions = '0 0 0'
    input_files = 'xolotl_userobj.i'
  [../]
[]
