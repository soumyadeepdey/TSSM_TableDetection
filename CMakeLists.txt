cmake_minimum_required(VERSION 2.8)
project( TableDetection )
set(BUILD_SHARED_LIBS OFF)

find_package( OpenCV REQUIRED )

include_directories(${OpenCV_INCLUDE_DIRS})

SET(TableDetection_HEADERS 	
				StandardHeaders.h
				sorting.h
				StatisticalFunctions.h
				folder.h
				binarization.h
				PixelValidation.h
				connectedcomponent.h				
				GeneralFunctions.h
				HSV.h
				Image_proc_functions.h
				Morphology.h				
				RectangleTest.h
				ScalarColorFeature.h
				SmoothingGapfilling.h								
				StrokeWidth.h
				TableDetection.h)  


SET(TableDetection_SOURCES  	
				TableDetection.cpp
				folder.cpp
				sorting.cpp
				StatisticalFunctions.cpp
				PixelValidation.cpp
				binarization.cpp
				connectedcomponent.cpp
				HSV.cpp
				GeneralFunctions.cpp				
				Image_proc_functions.cpp
				Morphology.cpp				
				RectangleTest.cpp
				ScalarColorFeature.cpp
				SmoothingGapfilling.cpp								
				StrokeWidth.cpp
				)
				


# Or set it globally for all targets:
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

ADD_EXECUTABLE(TabDet ${TableDetection_HEADERS} ${TableDetection_SOURCES})
TARGET_LINK_LIBRARIES(TabDet ${OpenCV_LIBS})
set(CMAKE_BUILD_TYPE Debug)

