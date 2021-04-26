# default length unit: nm
# default time unit: s
# default mass unit: ?
[Mesh]
  type = XolotlReflectedMesh
  dim = 2
  XolotlInput_path_name = './param_x2x_2D_noCnoR.txt'
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
  [./XolotlXeRate]
    order = FIRST
    family = LAGRANGE
  [../]
  [./XolotlXeMono]
    order = FIRST
    family = LAGRANGE
  [../]
  [./XolotlVolumeFraction]
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
    dt = 5000.0
  [../]
  start_time = 0
  end_time = 1.2e8
  num_steps = 10  
[]

[Outputs]
  exodus = true
[]

[MultiApps]
  [./sub_app]
    type = TransientMultiApp
    positions = '0 0 0'
input_files = 'x2x_subapp_2D_noCnoR.i'
app_type = coupling_xolotlApp
execute_on = TIMESTEP_END
library_path = 'lib'
[../]
[]

[Transfers]
[./fromsubrate]
type = MultiAppMeshFunctionTransfer
direction = from_multiapp
multi_app = sub_app
source_variable = Auxv
variable = XolotlXeRate
[../]
[./fromsubmono]
type = MultiAppMeshFunctionTransfer
direction = from_multiapp
multi_app = sub_app
source_variable = AuxMono
variable = XolotlXeMono
[../]
[./fromsubfrac]
type = MultiAppMeshFunctionTransfer
direction = from_multiapp
multi_app = sub_app
source_variable = AuxFrac
variable = XolotlVolumeFraction
[../]
[]
