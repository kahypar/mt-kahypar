function(hide_symbols target)
  if(NOT TARGET "${target}")
    message(FATAL_ERROR "hide_symbols: target '${target}' does not exist")
  endif()
  set_target_properties(${target} PROPERTIES
      C_VISIBILITY_PRESET hidden
      CXX_VISIBILITY_PRESET hidden
      VISIBILITY_INLINES_HIDDEN YES
  )
endfunction()

function(copy_dependencies target)
  if(NOT TARGET "${target}")
    message(FATAL_ERROR "copy_dependencies: target '${target}' does not exist")
  endif()

  if(WIN32)
    # copy (.dll) dependencies to target location
    add_custom_command(TARGET "${target}" POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              $<TARGET_RUNTIME_DLLS:${target}>
              $<TARGET_FILE_DIR:${target}>
      COMMAND_EXPAND_LISTS
    )
  endif()
endfunction()

function(add_mtkahypar_executable target)
  add_executable(${ARGV})
  get_target_property(is_alias ${target} ALIASED_TARGET)
  if(NOT is_alias)
    hide_symbols(${target})
  endif()
  copy_dependencies(${target})
endfunction()

function(add_mtkahypar_library target)
  add_library(${ARGV})
  get_target_property(is_alias ${target} ALIASED_TARGET)
  if(NOT is_alias)
    hide_symbols(${target})
  endif()
  copy_dependencies(${target})
endfunction()
