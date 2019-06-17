# default length unit: nm
# default time unit: s
# default mass unit: ?
[Mesh]
  type = PETScDMDAMesh
  dim = 2
  XolotlInput_path_name = './params_NE_2D_withPFsimple.txt'
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
[]

[Problem]
 type = XolotlProblem
 sync_variable = Auxv
 sync_GB = AuxGB
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
    dt = 5000.0
  [../]
  start_time = 0
  end_time = 7e7
  # end_time = 20000.0
[]
[Outputs]
  exodus = true
[]
