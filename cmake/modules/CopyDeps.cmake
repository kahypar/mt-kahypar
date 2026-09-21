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
