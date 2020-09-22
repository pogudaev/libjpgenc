project(jpgenc C)
cmake_minimum_required(VERSION 2.8)

option(BUILD_STATIC_LIBS "Build the static library" OFF)
option(BUILD_SHARED_LIBS "Build the shared library" ON)

#============================== VERSION ================================
set(LIB_VERSION_MAJOR 1)
set(LIB_VERSION_MINOR 0)
set(LIB_VERSION_PATCH 2)
set(LIB_VERSION_STRING ${LIB_VERSION_MAJOR}.${LIB_VERSION_MINOR}.${LIB_VERSION_PATCH})
#=======================================================================

#список хедеров для установки
set(CT_PUBLIC_HDRS 
    include/jpgenc.h
    #остальные хедеры не нужны для использования библиотеки
)

#указываем путь для хедеров
include_directories(BEFORE ./include)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -Werror -Wall -std=gnu99 -Os")
aux_source_directory(./src SRC_LIST)

if (BUILD_STATIC_LIBS OR BUILD_SHARED_LIBS)

	if (BUILD_STATIC_LIBS)
        message("build static lib")
		add_library(${PROJECT_NAME}_static STATIC ${SRC_LIST})
		install(TARGETS ${PROJECT_NAME}_static DESTINATION lib)
	endif (BUILD_STATIC_LIBS)

	if (BUILD_SHARED_LIBS)
        message("build shared lib")
		add_library(${PROJECT_NAME} SHARED ${SRC_LIST})
		install(TARGETS ${PROJECT_NAME} DESTINATION lib)
		set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${LIB_VERSION_STRING} SOVERSION ${LIB_VERSION_MAJOR})
	endif (BUILD_SHARED_LIBS)

	install(FILES ${CT_PUBLIC_HDRS} DESTINATION include/jpgenc)

else (BUILD_STATIC_LIBS OR BUILD_SHARED_LIBS)

	message(FATAL_ERROR "no target for build!")

endif (BUILD_STATIC_LIBS OR BUILD_SHARED_LIBS)

add_custom_target(uninstall
	COMMAND xargs rm -v < install_manifest.txt
)

