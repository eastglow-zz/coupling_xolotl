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
  library_path_name = '/Users/donguk.kim/projects/xolotl-build/lib/libxolotlInter.dylib'
  XolotlInput_path_name = './params_NE.txt'
[]
[Problem]
  kernel_coverage_check = false
  solve = false
[]
