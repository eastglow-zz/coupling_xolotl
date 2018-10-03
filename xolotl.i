[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
[]
[Variables]
  [./d]
  [../]
[]
[Executioner]
  type = XolotlExecutioner
  [./TimeStepper]
    type = XolotlTimeStepper
    # library_path_name ='/Users/sophie/GitWorkspace/xolotl-moose-source/build/lib/libxolotlInter.dylib'
    # library_path_name = '/Users/donguk.kim/projects/xolotl-build/lib/libxolotlInter.dylib'
  [../]
  end_time = 10
[]
[Problem]
  kernel_coverage_check = false
  solve = false
[]
