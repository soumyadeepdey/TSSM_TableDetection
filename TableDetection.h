

#include "StandardHeaders.h"
#include "sorting.h"
#include "StatisticalFunctions.h"
#include "PixelValidation.h"
#include "HSV.h"
#include "binarization.h"
#include "connectedcomponent.h"
#include "folder.h"
#include "GeneralFunctions.h"
#include "Image_proc_functions.h"
#include "Morphology.h"
#include "RectangleTest.h"
#include "ScalarColorFeature.h"
#include "SmoothingGapfilling.h"
#include "StandardHeaders.h"
#include "StrokeWidth.h"

typedef struct lineslopeintercept // y = mx + c ==> ay + bx = c ==> a = 1, b = -m, c
{
  double a;
  double b;
  double c;
}linesi;

namespace IITkgp_functions {
  
  /**
   * @function GetWordandLine
   * @param input 	vector<vector<Point> > contours = Set of obtained contours
   * 		 	vector<Rect> boundRect = Boundary rectangles obtained from the computed contours
   * 			vector<bool> validBlock =  Contours which are considered as valid for Table Detection
   * @brief		This function is used to compute lines from the valid connected components
   * @output 		vector<vector<int> > LineBlocknum ==> where LineBlocknum.size() = number of lines obtained from the given connected component information
   * 							      and  LineBlocknum[i].size() = Number of connected component present in that line
   */
  
vector<vector<int> >  GetWordandLine(vector<vector<Point> > contours,vector<Rect> boundRect,vector<bool> validBlock);



void sortcc(vector<int> & src, vector<int> & mask);




vector<vector<int> > SortccWithinLine(vector<vector<int> > LineBlocknum,vector<Rect> boundRect);


vector<vector<vector<RectDist> > > DistacceBetweenCCLine(vector<vector<int> > LineSortedBlocknum, vector<Rect> boundRect);


vector<vector<int> > JoinCCGap(vector<int> line, vector<vector<RectDist> > distance, int hgap, int vgap);


vector<vector<vector<int> > > GetJoinedCC(vector<vector<int> > line, vector<vector<vector<RectDist> > > distance, int hgap, int vgap);


vector<vector<int> > GetProbableTableLine(vector<vector<int> > mainline, vector<vector<vector<int> > > line);


  
}
