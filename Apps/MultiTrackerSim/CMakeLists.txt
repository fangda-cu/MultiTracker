append_files (Headers "h" . drivers drivers/curlnoise)
append_files (Sources "cpp" . drivers drivers/curlnoise)

include_directories (${EIGEN_INCLUDE_DIR})

add_executable (MultiTrackerSim ${Headers} ${Sources})
target_link_libraries (MultiTrackerSim LosTopos)

