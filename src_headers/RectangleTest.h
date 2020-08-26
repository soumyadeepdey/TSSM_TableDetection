


#include "StandardHeaders.h"


namespace IITkgp_functions {

  
    /**
     * @Function: PointRectangleTest
     * @brief : Take 1 rectangle and a Point as input
     * @brief : Test whether the Given Point is inside the Given Rectangle or Inside
     * @return : 	0: Point is Outside of Rectangle
     * 		1: Point is inside the given Rectangle
     * */

    int PointRectangleTest(Rect GivenRect, Point GivenPoint);

    /**
     * @Function: FindOverlappingRectangles
     * @brief : Take 2 rectangle as an input
     * @return : 	0: Rect 1 and Rect 2 are disjoint
     * 		1: Rect 1 is inside Rect 2
     * 	    	2: Rect 2 is inside Rect 1
     * 	    	3: Rect 1 and 2 are partly overlapped
     *
     *
     * */

    int FindOverlappingRectangles(Rect first, Rect second);

  
}