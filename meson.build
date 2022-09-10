project('olcPGEX_ViewPort', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3'])

external_deps = [dependency('threads'), dependency('gl')]
cpp = meson.get_compiler('cpp')

if host_machine.system() == 'windows'
  # no platform specific deps are needed
elif host_machine.system() == 'darwin'
  external_deps += [dependency('libpng'),
                    dependency('appleframeworks', modules: ['GLUT'])]
else
  external_deps += [dependency('libpng'),
                    dependency('x11'),
                    cpp.find_library('stdc++fs', required: false)]
endif

exe = executable('olcpgex_viewport', 'example.cpp',
  install : true,
  dependencies : external_deps,
)

test('basic', exe)
