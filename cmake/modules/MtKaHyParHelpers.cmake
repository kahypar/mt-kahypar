function(hide_symbols target)
  if(NOT TARGET "${target}")
    message(FATAL_ERROR "hide_symbols: target '${target}' does not exist")
  endif()

  get_target_property(is_alias ${target} ALIASED_TARGET)
  if(NOT is_alias)
    set_target_properties(${target} PROPERTIES
        C_VISIBILITY_PRESET hidden
        CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN YES
    )
  endif()
endfunction()

function(copy_dependencies target)
  if(NOT TARGET "${target}")
    message(FATAL_ERROR "copy_dependencies: target '${target}' does not exist")
  endif()

  get_target_property(is_alias ${target} ALIASED_TARGET)
  get_target_property(target_type ${target} TYPE)
  if(WIN32 AND (NOT is_alias) AND
      (target_type MATCHES "^(EXECUTABLE|(STATIC|SHARED|MODULE)_LIBRARY)$"))
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
  hide_symbols(${target})
  copy_dependencies(${target})
endfunction()

function(add_mtkahypar_library target)
  add_library(${ARGV})
  hide_symbols(${target})
  copy_dependencies(${target})
endfunction()
