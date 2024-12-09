# Check if the required APIs for thread pinning are available

set(THREAD_PINNING_WORKS FALSE)
if(UNIX)
	# Check if sched_getcpu/sched_setaffinity are available
	check_include_files("sched.h" HEADER_AVAILABLE LANGUAGE CXX)
	if(${HEADER_AVAILABLE})
		# we need to actually check the available symbols
		include(CheckSourceCompiles)
		check_source_compiles(CXX
		"#include <sched.h>
		int main() {
			int cpu_id = sched_getcpu();
			int num_cpus = 4;
			const size_t size = CPU_ALLOC_SIZE(num_cpus);
			cpu_set_t mask;
			CPU_ZERO(&mask);
			CPU_SET(cpu_id, &mask);
			int err = sched_setaffinity(0, size, &mask);
			return 0;
		}"
		EXAMPLE_COMPILES)
		if(${EXAMPLE_COMPILES})
			set(THREAD_PINNING_WORKS TRUE)
		endif()
	endif()
elseif(WIN32)
	set(THREAD_PINNING_WORKS TRUE)
endif()
