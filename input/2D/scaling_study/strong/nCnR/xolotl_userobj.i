# default length unit: nm
# default time unit: s
# default mass unit: ?
[Mesh]
  type = XolotlReflectedMesh
  dim = 2
  XolotlInput_path_name = './params_NE_2D_20umsq_401sq.txt'
[]

[AuxVariables]
  [./Auxv]
    order = FIRST
    family = LAGRANGE
  [../]
  [./AuxGB]
    order = FIRST
    family = LAGRANGE
  [../]
  [./AuxMono]
    order = FIRST
    family = LAGRANGE
  [../]
  [./AuxFrac]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Problem]
 type = XolotlProblem
 sync_rate = Auxv
 sync_GB = AuxGB
 sync_mono = AuxMono
 sync_frac = AuxFrac
[]

[Variables]
  [./d]
  [../]
[]
[ICs]
  [./Init_Aux_const]
    type = ConstantIC
    variable = Auxv
    value = 0
  [../]
  [./Init_Aux_gb]
    type = ConstantIC
    variable = AuxGB
    value = 1
  [../]
[]

[Executioner]
  type = Transient
  [./TimeStepper]
    type = ConstantDT
    dt = 0.5
  [../]
  start_time = 0
  end_time = 1e8
  # end_time = 20000.0
[]
[Outputs]
  exodus = true
[]
