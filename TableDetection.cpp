



#include "TableDetection.h"

RNG rng(12345);

using namespace IITkgp_functions;


vector<vector<int> >  IITkgp_functions::GetWordandLine(vector<vector<Point> > contours,vector<Rect> boundRect,vector<bool> validBlock)
{
  vector<vector<int> > LineFormed;
  
  for(int i =0;i<contours.size();i++)
  {    
    if(validBlock[i])
    {
      vector<int> tempw;
    
      Rect TempRect;
      TempRect = boundRect[i];
      Point2i tc1 = TempRect.tl();
      Point2i bc1 = TempRect.br();
      int miny,maxy;
      miny = tc1.y;
      maxy = bc1.y;
      for(int j=0;j<contours.size();j++)
      {
	if(validBlock[j] )
	{
	  Rect TempRect1;
	  TempRect1 = boundRect[j];
	  Point2i tc2 = TempRect1.tl();
	  Point2i bc2 = TempRect1.br();	
	  if(tc2.y >= miny && tc2.x <= maxy )// if 2nd top left point within miny and maxy
	  {
	    tempw.push_back(j);

	    if(miny > tc2.y)
	    miny = tc2.y;
	    if(maxy < bc2.y)
	      maxy = bc2.y;
	  }
	  else if(bc2.y >= miny && bc2.y <= maxy)// if 2nd bottom right point within miny and maxy
	  {
	    tempw.push_back(j);

	    if(miny > tc2.y)
	    miny = tc2.y;
	    if(maxy < bc2.y)
	      maxy = bc2.y;
	  }
	  else if(miny >= tc2.x && maxy <=bc2.x)
	  {
	    tempw.push_back(j);

	    if(miny > tc2.y)
	    miny = tc2.y;
	    if(maxy < bc2.y)
	      maxy = bc2.y;
	  }
	}
      }
      LineFormed.push_back(tempw);
    }
  }
  
  return(LineFormed);
  
}





void IITkgp_functions::sortcc(vector<int> & src, vector<int> & mask)
{
  int length = mask.size();

    for (int i = 0; i < length; ++i)
    {
        bool swapped = false;
        for (int j = 0; j < length - (i+1); ++j)
        {
            if (mask[j] > mask[j+1])
            {
                swaping<int>(mask, j, j+1);
		swaping<int>(src, j, j+1);
                swapped = true;
            }
        }
        
        if (!swapped) break;
    }
}

void sortCLines(vector<Vec4i> & src, vector<int> & mask)
{
  int length = mask.size();

    for (int i = 0; i < length; ++i)
    {
        bool swapped = false;
        for (int j = 0; j < length - (i+1); ++j)
        {
            if (mask[j] > mask[j+1])
            {
                swaping<int>(mask, j, j+1);
		Vec4i tmp = src[j];
		src[j] = src[j+1];
		src[j+1] = tmp;
                swapped = true;
            }
        }
        
        if (!swapped) break;
    }
}


vector<vector<int> > IITkgp_functions::SortccWithinLine(vector<vector<int> > LineBlocknum,vector<Rect> boundRect)
{
  for(int i=0;i<LineBlocknum.size();i++)
 {
   vector<int> LineCC = LineBlocknum[i];
   vector<int> CCXCord;
   
   for(int j=0;j<LineCC.size();j++)
   {
     Rect temprect = boundRect[LineCC[j]];
     Point2i tl = temprect.tl();
     CCXCord.push_back(tl.x);
   }
   sortcc(LineCC,CCXCord);
   LineBlocknum[i] = LineCC;
   LineCC.clear();
   CCXCord.clear();
 }
 return(LineBlocknum);
}









vector<vector<vector<RectDist> > > IITkgp_functions::DistacceBetweenCCLine(vector<vector<int> > LineSortedBlocknum, vector<Rect> boundRect)
{
  vector<vector<vector<RectDist> > > CCDistanceLine;
  for(int i=0;i<LineSortedBlocknum.size();i++)
  {
    vector<int> Sortedcc = LineSortedBlocknum[i];
    vector<vector<RectDist> > Distance;
    for(int j=0;j<Sortedcc.size();j++)
    {
      vector<RectDist> dist;
      for(int k=0;k<Sortedcc.size();k++)
      {
	dist.push_back(DistanceBetweenRectangle(boundRect[Sortedcc[j]],boundRect[Sortedcc[k]]));
      }
      Distance.push_back(dist);
      dist.clear();
     //Distance.push_back(DistanceBetweenRectangle(boundRect[Sortedcc[j]],boundRect[Sortedcc[j+1]])); 
    }
    CCDistanceLine.push_back(Distance);
    Distance.clear();
    Sortedcc.clear();
  }
  return(CCDistanceLine);
}







vector<vector<int> > IITkgp_functions::JoinCCGap(vector<int> line, vector<vector<RectDist> > distance, int hgap, int vgap)
{
  vector<vector<int> > joindedcc;  
  vector<bool> flag(line.size(),true);
  
  for(int i=0;i<line.size();i++)
  {
    if(flag[i])
    {
      vector<int> lineWords;
      lineWords.push_back(i);
      flag[i] = false;
      for(int j=0;j<line.size();j++)
      {
	if(flag[j])
	{ 
	  if(distance[i][j].vdist < vgap)
	  {
	    if(distance[i][j].hdist < hgap)
	    {
	      lineWords.push_back(j);
	      flag[j] = false;
	    }  
	  }
	  else
	  {
	    if(distance[i][j].hdist < hgap)
	    {
	      lineWords.push_back(j);
	      flag[j] = false;
	    }
	  }
        }
      }
      joindedcc.push_back(lineWords);
    }
  }
  
  
  return(joindedcc);
}


 vector<vector<vector<int> > > IITkgp_functions::GetJoinedCC(vector<vector<int> > line, vector<vector<vector<RectDist> > > distance, int hgap, int vgap)
{
  vector<vector<vector<int> > > JoinedLineBlock;
 
  for(int i=0;i<line.size();i++)
  {
    vector<int> templine = line[i];
    vector<vector<RectDist> > dis = distance[i];
    vector<vector<int> > tempjl = JoinCCGap(templine,dis,hgap,vgap);
    dis.clear();
    templine.clear();
    JoinedLineBlock.push_back(tempjl);
    tempjl.clear();
  }  
  return(JoinedLineBlock);
}










vector<vector<int> > IITkgp_functions::GetProbableTableLine(vector<vector<int> > mainline, vector<vector<vector<int> > > line)
{
  vector<vector<int> > Probtabline;
  
  for(int i=0;i<line.size();i++)
  {
    vector<vector<int> > tmpline = line[i];
    if(!tmpline.empty())
    {
      if(tmpline.size()>1)
      {
	Probtabline.push_back(mainline[i]);
      }
    }
  }
  
  return(Probtabline);
}

int NewPointRectangleTest(Rect GivenRect, Point GivenPoint, float th)
{
  Point tl,br;
  tl = GivenRect.tl();
  br = GivenRect.br();
  int flag;
  /*
  if((GivenPoint.x>=tl.x && GivenPoint.x<=br.x) && (GivenPoint.y<=tl.y && GivenPoint.y>=br.y))
  {
    flag = 1;
    printf("point inside\n");
    return(flag);
  }
  */
  if((GivenPoint.x>=tl.x && GivenPoint.x<=br.x) && (GivenPoint.y>=tl.y && GivenPoint.y<=br.y))
  {
    flag = 1;
    //printf("point inside\n");
    return(flag);
  }
  else if(GivenPoint.x>=tl.x && GivenPoint.x<=br.x && (abs(GivenPoint.y-tl.y)*1.0 < th || abs(GivenPoint.y-br.y)*1.0 < th))
  {
    flag = 1;
    //printf("point inside\n");
    return(flag);
  }
  else if(GivenPoint.y>=tl.y && GivenPoint.y<=br.y && (abs(GivenPoint.x-tl.x)*1.0 < th || abs(GivenPoint.x-br.x)*1.0 < th))
  {
    flag = 1;
    //printf("point inside\n");
    return(flag);
  }
  else
  {
    flag = 0;
    return(flag);
  }
}


double CalculateSlopeofaLine(Vec4i line)
{
  double y = (line[3]-line[1]);
  double x = (line[2] -line[0]);
  double slope = y/x;
  return(slope);
}

int CalculateyIntercept(Vec4i line)
{
  double m = CalculateSlopeofaLine(line);
  double y = line[3]*1.0;
  double x = line[2]*1.0;
  double c = y - (m*x);
  double temp = c + 0.5;
  int yintercept =(int) floor(temp);
  return(yintercept);
}

int CalculatexIntercept(Vec4i line)
{
  double m = CalculateSlopeofaLine(line);
  double y = line[3]*1.0;
  double x = line[2]*1.0;
  double c = y - (m*x);
  double a = (c/m)*(-1.0);
  double temp = a + 0.5;
  int xintercept = (int) floor(a);
  return(xintercept);
}

double FindAngleofaLine(Vec4i line)
{    
    double y = (line[3]-line[1]);
    double x = (line[2] -line[0]);
    double angle;
    angle = atan2(y,x);
    angle = (angle*360)/(2*PI);
    return(angle);
}


/**
 * return 1 if horizontal
 * return 2 if vertical
 */
int LineType(Vec4i line)
{
    int type = 0;
    double y = (line[3]-line[1]);
    double x = (line[2] -line[0]);
    double angle;
    angle = atan2(y,x);
    angle = (angle*360)/(2*PI);
    angle = abs(angle);
    if(angle < 45) 
      type = 1; //Horizontal LINE
    if(angle >=45)
      type = 2; // Vertical
    
    return(type);
}


double DistanceBetweentwoParallelLine(Vec4i line1, Vec4i line2)
{
  
  
  
  double Dist;
   
   double slope1 = CalculateSlopeofaLine(line1);
   if(line1[0] == line1[2])
  {
    int xintercept1 = line1[0];
    int xintercept2 = line2[0];
    Dist = abs(xintercept1 - xintercept2);
  }
  else if(line1[1] == line1[3])
  {
    int yintercept1 = line1[1];
    int yintercept2 = line2[1];
    Dist = abs(yintercept1 - yintercept2);
  }
  else
  {
    int yintercept1 = CalculateyIntercept(line1);
    int yintercept2 = CalculateyIntercept(line2);
    /*
    double a1 = 1;
    double b1 = slope1 * -1;;
    double c1 = yintercept1 * -1;
    double c2 = yintercept2 * -1;
    Dist = abs(c2 - c1);
    double temp = (a1*a1)+(b1*b1);
    temp = sqrt(temp);
    Dist = Dist/temp;
    */
    Dist = abs(yintercept2-yintercept1);
    Dist = Dist/sqrt((slope1*slope1)+1);
  }
 // printf("Distance = %lf\n",Dist);
  
 /* if(LineType(line1)==1)
  {
    int xintercept1 = line1[0];
    int xintercept2 = line2[0];
    Dist = abs(xintercept1 - xintercept2);
  }
  else
  {
    int yintercept1 = line1[1];
    int yintercept2 = line2[1];
    Dist = abs(yintercept1 - yintercept2);
  }
  */
  return(Dist);
}

vector<vector<Vec4i> > FindParallelLines(vector<Vec4i> lines,double threshold)
{
  vector<double> slopes;
  
  for(int i=0;i<lines.size();i++)
  {
    double angle;
    angle = FindAngleofaLine(lines[i]);
    slopes.push_back(angle);
  }
  
  vector<vector<Vec4i> > ParallelLineSet;
  vector<bool> flag(lines.size(),true);
  
  for(int i=0;i<lines.size();i++)
  {
    if(flag[i])
    {
      flag[i] = false;
      vector<Vec4i> ParallelLine;
      ParallelLine.push_back(lines[i]);
      for(int j=0;j<lines.size();j++)
      {
	if(flag[j])
	{
	  if(abs(slopes[i]-slopes[j]) <= threshold)
	  {
	    flag[j] = false;
	    ParallelLine.push_back(lines[j]);
	  }
	  else if(abs(slopes[i]-slopes[j]) > 90)
	  {
	    if(180 - abs(slopes[i]-slopes[j]) <= threshold)
	    {
	      flag[j] = false;
	      ParallelLine.push_back(lines[j]);
	    }
	  }
	}
      }
      if(!ParallelLine.empty())
      {
	ParallelLineSet.push_back(ParallelLine);
      }
    }
  }
  //cout << "parallel Line set\t" << ParallelLineSet.size() << endl;
  return(ParallelLineSet);
  
}

bool isEqual(const Vec4i& _l1, const Vec4i& _l2)
{
    Vec4i l1(_l1), l2(_l2);

    float length1 =(float) sqrt((l1[2] - l1[0])*(l1[2] - l1[0]) + (l1[3] - l1[1])*(l1[3] - l1[1]));
    float length2 =(float) sqrt((l2[2] - l2[0])*(l2[2] - l2[0]) + (l2[3] - l2[1])*(l2[3] - l2[1]));

    float product = (l1[2] - l1[0])*(l2[2] - l2[0]) + (l1[3] - l1[1])*(l2[3] - l2[1]);

    if (fabs(product / (length1 * length2)) < cos(CV_PI / 30))
        return false;

    float mx1 = (l1[0] + l1[2]) * 0.5f;
    float mx2 = (l2[0] + l2[2]) * 0.5f;

    float my1 = (l1[1] + l1[3]) * 0.5f;
    float my2 = (l2[1] + l2[3]) * 0.5f;
    float dist =(float) sqrt((mx1 - mx2)*(mx1 - mx2) + (my1 - my2)*(my1 - my2));
   // cout << dist << endl;

    if (dist > std::max(length1, length2) * 0.5f)
        return false;

    return true;
}


vector<vector<vector<Vec4i> > > FindConnectedLines(vector<vector<Vec4i> > ParallelLineSet)
{
  vector<vector<vector<Vec4i> > > ConnectedParallelLines;
  for(int i=0;i<ParallelLineSet.size();i++)
  {
    vector<Vec4i> lines = ParallelLineSet[i];
    vector<int> labels;
    int numberOfLines = partition(lines, labels, isEqual);
    vector<vector<Vec4i> > Clines(numberOfLines);
    for(int j=0;j<labels.size();j++)
    {
      Clines[labels[j]].push_back(lines[j]);
    }
    ConnectedParallelLines.push_back(Clines);
  }
  return ConnectedParallelLines;
}

float lineEqualityThreshold;

bool isEqualNew(const Vec4i& _l1, const Vec4i& _l2)
{
    Vec4i l1(_l1), l2(_l2);

    float dist;
    
    if(LineType(l1)==2)// vertical
    {
	float midx1 = (l1[0] + l1[2])*0.5;
	float midx2 = (l2[0]+l2[2])*0.5;
	dist =(float) abs(midx1 - midx2); 
    }
    else if(LineType(l2)==1)//horizontal
    {
	float midy1 = (l1[1] + l1[3])*0.5;
	float midy2 = (l2[1]+l2[3])*0.5;
	dist =(float) abs(midy1 - midy2); 
	//cout << dist << endl;
    }
    //cout << dist << endl;
    if(dist > lineEqualityThreshold)
      return false;
    else
      return true;
}


vector<vector<vector<Vec4i> > > FindConnectedLinesNew(vector<vector<Vec4i> > ParallelLineSet,double threshold)
{
  vector<vector<vector<Vec4i> > > ConnectedParallelLines;
  lineEqualityThreshold =(float) threshold;
  for(int i=0;i<ParallelLineSet.size();i++)
  {
    vector<Vec4i> lines = ParallelLineSet[i];
    vector<int> labels;
    int numberOfLines = partition(lines, labels, isEqualNew);
    //cout << "numberOfLines=" << numberOfLines << endl;
    vector<vector<Vec4i> > Clines(numberOfLines);
    for(int j=0;j<labels.size();j++)
    {
      Clines[labels[j]].push_back(lines[j]);
    }
    ConnectedParallelLines.push_back(Clines);
  }
  return ConnectedParallelLines;
}



vector<vector<vector<Vec4i> > > SplitConnectedLines(vector<vector<vector<Vec4i> > > ConnectedParallelLines, float threshold)
{
  lineEqualityThreshold = threshold;
  vector<vector<vector<Vec4i> > > ConnectedLines;
  for(int i=0;i<ConnectedParallelLines.size();i++)
  {
    vector<vector<Vec4i> > TempLines;
    for(int j=0;j<ConnectedParallelLines[i].size();j++)
    {
      vector<Vec4i> lines = ConnectedParallelLines[i][j];
      int ltype = LineType(lines[0]);
      vector<int> mask;
      vector<int> dist;
      if(ltype == 1)
      {
	//cout << "Horizontal Line Pline=" << i << "Cline" << j << endl;
	for(int k=0;k<lines.size();k++)
	{
	  Vec4i line = lines[k]; 
	  mask.push_back(min(line[0],line[2]));
	}
	sortCLines(lines,mask);
	dist.push_back(0);
	int max_xb = max(lines[0][0],lines[0][2]);
	for(int k=1;k<lines.size();k++)
	{
	  Vec4i line2 = lines[k];
	  Vec4i line1 = lines[k-1];
	  max_xb = max(max_xb,max(line1[2],line1[0]));
	  dist.push_back( min(line2[0],line2[2]) - max_xb);
	  //cout << dist[k] << endl;
	  //dist.push_back(min((line2[0] - line1[2]),(line2[0] - line1[0])));
	}
      }
      else if(ltype == 2)
      {
	//cout << "Vertical Line Pline=" << i << "Cline" << j << endl;
	for(int k=0;k<lines.size();k++)
	{
	  Vec4i line = lines[k]; 
	  mask.push_back(min(line[1],line[3]));
	  //if(line[1] > line[3])
	   // cout << "My assumption is wrong" << endl;
	}
	sortCLines(lines,mask);
	dist.push_back(0);
	int max_yb = max(lines[0][3],lines[0][1]);
	for(int k=1;k<lines.size();k++)
	{
	  Vec4i line2 = lines[k];
	  Vec4i line1 = lines[k-1];
	 // int tv = max(line1[3],line1[1]);
	  max_yb = max(max_yb,max(line1[3],line1[1]));
	  //tv = min(line2[1],line2[3]);
	  dist.push_back( min(line2[1],line2[3]) - max_yb);
	 // cout << dist[k] << endl;
	 // dist.push_back(line2[1] - max_yb);
	}
      }
      mask.clear();
      vector<int> labels(dist.size());
      labels[0] = 0;
      
      for(int k=1;k<dist.size();k++)
      {
	if(dist[k] < threshold)
	  labels[k] = labels[k-1];
	else
	  labels[k] = labels[k-1] + 1;
      }
      vector<vector<Vec4i> > Clines(labels[labels.size()-1]+1);
      //cout << "Number of Lines = " << Clines.size() << endl;
      for(int k=0;k<labels.size();k++)
      {
	Clines[labels[k]].push_back(lines[k]);
      }
      for(int k=0;k<Clines.size();k++)
	TempLines.push_back(Clines[k]);
      Clines.clear();
    }
    ConnectedLines.push_back(TempLines);
  }
  
  return ConnectedLines;
}

vector<vector<vector<Vec4i> > > FindConnectedLines(vector<vector<Vec4i> > ParallelLineSet,double threshold)
{
 
  
  vector<vector<vector<Vec4i> > > ConnectedParallelLines(ParallelLineSet.size());
  
  for(int i=0;i<ParallelLineSet.size();i++)
  {
    
    vector<bool> flag(ParallelLineSet[i].size(),true);
    for(size_t j=0;j<ParallelLineSet[i].size();j++)
    {
      if(flag[j])
      {
	vector<Vec4i> tempclines;
	tempclines.push_back(ParallelLineSet[i][j]);
	flag[j] = false;
	for(size_t k=0;k<ParallelLineSet[i].size();k++)
	{
	  if(flag[k])
	  {
	   // printf("Distance Between %d\t%d\t%d\t",i,j,k);
	    double dist = DistanceBetweentwoParallelLine(ParallelLineSet[i][j],ParallelLineSet[i][k]);
	    if(dist <= threshold)
	    {
	      tempclines.push_back(ParallelLineSet[i][k]);
	      flag[k] = false;
	    }
	  }
	}
	ConnectedParallelLines[i].push_back(tempclines);
	tempclines.clear();	
      }     
    }
    flag.clear();    
  }
  
  return(ConnectedParallelLines);
  
}



vector<vector<Vec4i> > JoinedConnectedLineNew(vector<vector<vector<Vec4i> > > ConnectedParallelLines)
{
  vector<vector<Vec4i> > JoinedLines(ConnectedParallelLines.size());
  for( size_t k = 0; k < ConnectedParallelLines.size(); k++ )
  {
    vector<vector<Vec4i> > ConnectedLines = ConnectedParallelLines[k];
    for(size_t i = 0; i < ConnectedLines.size(); i++)
    {
      vector<Vec4i> lines = ConnectedLines[i];
      Vec4i line;
      double angle;
      angle = abs(FindAngleofaLine(lines[0]));
      if(LineType(lines[0])==1)//horizontal
      {
	int xl,xr,y;
	xl = min(lines[0][0],lines[0][2]); xr = max(lines[0][0],lines[0][2]); y = min(lines[0][1],lines[0][3]);
	for(int j=0;j<lines.size();j++)
	{
	    xl = min(min(lines[j][0],lines[j][2]),xl);
	    xr = max(max(lines[j][0],lines[j][2]),xr);
	    y = min(min(lines[j][1],lines[j][3]),y);
	    //cout << xl << "," << y << ": " << xr << ","  << y << endl;
	}
	line[0] = xl; line[2] = xr; line[1] = y; line[3] = y;
	//if(angle < 1.5)
	  JoinedLines[k].push_back(line);
      }
      else if(LineType(lines[0])==2)
      {
	int yt,yb,x;
	yt = min(lines[0][1],lines[0][3]); yb = max(lines[0][1],lines[0][3]); x = min(lines[0][0],lines[0][2]);
	
	for(int j=0;j<lines.size();j++)
	{
	    yt = min(min(lines[j][1],lines[j][3]),yt);
	    yb = max(max(lines[j][1],lines[j][3]),yb);
	    x = min(min(lines[j][0],lines[j][2]),x);
	   // cout << x << "," << yt << ": " << x << ","  << yb << endl;
	}
	line[0] = x; line[2] = x; line[1] = yt; line[3] = yb;
	//cout << "vertical line length is " << yb-yt << endl;
	//if((90 - angle) < 1.5 )
	  JoinedLines[k].push_back(line);
      }
      
    }
  }
  
  return JoinedLines;
}


vector<vector<Vec4i> > JoinedConnectedLine(vector<vector<vector<Vec4i> > > ConnectedParallelLines)
{
  vector<vector<Vec4i> > JoinedLines(ConnectedParallelLines.size());
  
  for( size_t k = 0; k < ConnectedParallelLines.size(); k++ )
  {
    vector<vector<Vec4i> > ConnectedLines = ConnectedParallelLines[k];
    for(size_t i = 0; i < ConnectedLines.size(); i++)
    {
      Vec4i JoinedLine;
      vector<Vec4i> TLines = ConnectedLines[i];
      int type = LineType(TLines[0]);
      double angle;
      angle = abs(FindAngleofaLine(TLines[0]));
      Point2i P,Q;
      P.x = min(TLines[0][0],TLines[0][2]);
      if(P.x == TLines[0][0])
      {
	P.y = TLines[0][1];
	Q.x = TLines[0][2];
	Q.y = TLines[0][3];
      }
      else
      {
	P.y = TLines[0][3];
	Q.x = TLines[0][0];
	Q.y = TLines[0][1];
      }
	
      if(type == 1) // Horizontal line
      {
	if(angle < 0.5)
	{
	  for(size_t j = 0; j < TLines.size(); j++)
	  {
	    Point TP,TQ;
	    TP.x = min(TLines[j][0],TLines[j][2]);
	    if(TP.x == TLines[j][0])
	    {
	      TP.y = TLines[j][1];
	      TQ.x = TLines[j][2];
	      TQ.y = TLines[j][3];
	    }
	    else
	    {
	      TP.y = TLines[j][3];
	      TQ.x = TLines[j][0];
	      TQ.y = TLines[j][1];
	    }
	    
	    P.x = min(TP.x,P.x);
	    
	    if(P.x == TP.x)
	    {
	      P.y = TP.y;
	    }
	    
	    Q.x = max(TQ.x,Q.x);
	    if(Q.x == TQ.x)
	      Q.y = TQ.y;
	    
	  }
	  Q.y = P.y;
	  JoinedLine[0] = P.x; JoinedLine[1] = P.y;
	  JoinedLine[2] = Q.x; JoinedLine[3] = Q.y;
	  
	  
	  JoinedLines[k].push_back(JoinedLine);
	}
	else
	{
	  /*
	  double AvgAngle = 0.0;
	  for(size_t j = 0; j < TLines.size(); j++)
	  {
	    AvgAngle = AvgAngle + abs(FindAngleofaLine(TLines[j]));
	  }
	  AvgAngle = AvgAngle/TLines.size();
	  double m = tan(AvgAngle);
	  */
	  for(size_t j = 0; j < TLines.size(); j++)
	  {
	    Point TP,TQ;
	    TP.x = min(TLines[j][0],TLines[j][2]);
	    if(TP.x == TLines[j][0])
	    {
	      TP.y = TLines[j][1];
	      TQ.x = TLines[j][2];
	      TQ.y = TLines[j][3];
	    }
	    else
	    {
	      TP.y = TLines[j][3];
	      TQ.x = TLines[j][0];
	      TQ.y = TLines[j][1];
	    }
	    
	    P.x = min(TP.x,P.x);
	    
	    if(P.x == TP.x)
	    {
	      P.y = TP.y;
	    }
	    
	    Q.x = max(TQ.x,Q.x);
	    if(Q.x == TQ.x)
	      Q.y = TQ.y;
	  }
	  
	  double m = CalculateSlopeofaLine(TLines[0]);
	  int c = CalculateyIntercept(TLines[0]);
	  double ty = (m*Q.x)+c;
	  ty = ty + 0.5;
	  Q.y = (int) floor(ty);
	  
	  JoinedLine[0] = P.x; JoinedLine[1] = P.y;
	  JoinedLine[2] = Q.x; JoinedLine[3] = Q.y;
	  
	  JoinedLines[k].push_back(JoinedLine);
	}
      }
      else if(type == 2) // Vertical Line
      {
	if((90 - angle) < 0.5 )
	{
	  for(size_t j = 0; j < TLines.size(); j++)
	  {
	    Point TP,TQ;
	    TP.y = min(TLines[j][1],TLines[j][3]);
	    if(TP.y == TLines[j][1])
	    {
	      TP.x = TLines[j][0];
	      TQ.x = TLines[j][2];
	      TQ.y = TLines[j][3];
	    }
	    else
	    {
	      TP.x = TLines[j][2];
	      TQ.x = TLines[j][0];
	      TQ.y = TLines[j][1];
	    }
	    
	    P.y = min(TP.y,P.y);
	    
	    if(P.y == TP.y)
	    {
	      P.x = TP.x;
	    }
	    
	    Q.y = max(TQ.y,Q.y);
	    if(Q.y == TQ.y)
	      Q.x = TQ.x;
	    
	  }
	  Q.x = P.x;
	  JoinedLine[0] = P.x; JoinedLine[1] = P.y;
	  JoinedLine[2] = Q.x; JoinedLine[3] = Q.y;
	  
	  
	  JoinedLines[k].push_back(JoinedLine);
	}
	else
	{
	  /*
	  double AvgAngle = 0.0;
	  for(size_t j = 0; j < TLines.size(); j++)
	  {
	    AvgAngle = AvgAngle + abs(FindAngleofaLine(TLines[j]));
	  }
	  AvgAngle = AvgAngle/TLines.size();
	  double m = tan(AvgAngle);
	  */
	  for(size_t j = 1; j < TLines.size(); j++)
	  {
	    Point TP,TQ;
	    TP.y = min(TLines[j][1],TLines[j][3]);
	    if(TP.y == TLines[j][1])
	    {
	      TP.x = TLines[j][0];
	      TQ.x = TLines[j][2];
	      TQ.y = TLines[j][3];
	    }
	    else
	    {
	      TP.x = TLines[j][2];
	      TQ.x = TLines[j][0];
	      TQ.y = TLines[j][1];
	    }
	    
	    P.y = min(TP.y,P.y);
	    
	    if(P.y == TP.y)
	    {
	      P.x = TP.x;
	    }
	    
	    Q.y = max(TQ.y,Q.y);
	    if(Q.y == TQ.y)
	      Q.x = TQ.x;
	    
	  }
	  
	  double m = CalculateSlopeofaLine(TLines[0]);
	  int c = CalculateyIntercept(TLines[0]);
	 // double tx = P.x - (((P.y-Q.y)*1.0)/m);
	  double tx = ((Q.y - c)*1.0)/m;
	  tx = tx + 0.5;
	  Q.x = (int) floor(tx);
	  
	  JoinedLine[0] = P.x; JoinedLine[1] = P.y;
	  JoinedLine[2] = Q.x; JoinedLine[3] = Q.y;
	  
	  JoinedLines[k].push_back(JoinedLine);
	}
      }
      else
      {
	printf("ERROR: TYPE OF LINE IS INCORRECT\n");
	exit(0);
      }
    }
  }
  
  return(JoinedLines);
  
}

/*
vector<Vec4i> JoinedConnectedLine(vector<vector<Vec4i> > ConnectedLines)
{
  vector<Vec4i> JoinedLines;
  
  for( size_t i = 0; i < ConnectedLines.size(); i++ )
  {
    vector<Vec4i> TLines = ConnectedLines[i];
    Vec4i jl;
    jl[0] = min(TLines[0][0],TLines[0][2]);
    jl[2] = max(TLines[0][0],TLines[0][2]);
    jl[1] = min(TLines[0][1],TLines[0][3]);
    jl[3] = max(TLines[0][1],TLines[0][3]);   
    for(size_t j = 0; j < TLines.size(); j++)
    {
      int t;
      t=min(TLines[i][0],TLines[i][2]);
      if(jl[0] > t)
	jl[0] = t;
      t = max(TLines[i][0],TLines[i][2]);
      if(jl[2] < t)
	jl[2] = t;
      t = min(TLines[i][1],TLines[i][3]);
      if(jl[1] > t)
	jl[1] = t;
      t = max(TLines[i][1],TLines[i][3]);
      if(jl[3] < t)
	jl[3] = t;
    }
    double angle;
    angle = abs(FindAngleofaLine(TLines[0]));
    if(angle == 90)
    {
      jl[2] = jl[0];
    }
    else if(angle == 0)
    {
      jl[3] = jl[1];
    }
    else
    {
      int type = LineType(TLines[0]);
      double m = ((TLines[0][3] - TLines[0][1])*1.0)/((TLines[0][2] - TLines[0][0])*1.0);
      double c = CalculateyIntercept(TLines[0]);
      if(type ==  1) // Line type Horizontal
      {
	double tt = m*jl[2]+c;
	tt = tt + 0.5;
	jl[3] = (int) floor(tt);
      }
      else if(type ==  2) // Line typr Vertical
      {
	double tt = (jl[3]-c)/m;
	tt = tt + 0.5;
	jl[2] = (int) floor(tt);
      }
    }
    JoinedLines.push_back(jl);
    TLines.clear();
  }
  
  return(JoinedLines);
  
}

*/

linesi GetLinefrom2Point(Vec4i line1)
{
  double m = CalculateSlopeofaLine(line1);
  double c = CalculateyIntercept(line1);
  linesi T;
  T.a = 1.0;
  T.b = (m * -1.0);
  T.c = c;
  return(T);
}

double DistanceBeteenTwoPoint(Point2i P, Point2i Q)
{
  double dist;
  dist = (((P.x - Q.x)*(P.x - Q.x))+((P.y - Q.y)*(P.y - Q.y)))*1.0;
  dist = sqrt(dist);
  return(dist);
}

linesi GetPerpendicularLinePassingThroughaPoint(linesi L, Point2i P)
{
  linesi nl;
  nl.a = 1;
  nl.b = (1/L.b)*(-1.0);
  nl.c = (nl.a*P.y) + (nl.b*P.x);
  return(nl);
}

Point2i GetIntersectingPoint(linesi line1, linesi line2)
{
  Mat A = Mat(2,2,CV_64FC1);
  Mat C = Mat(2,1,CV_64FC1);
  
  A.at<double>(0,0) = line1.a;
  A.at<double>(0,1) = line1.b;
  C.at<double>(0,0) = line1.c;
  A.at<double>(1,0) = line2.a;
  A.at<double>(1,1) = line2.b;
  C.at<double>(1,0) = line2.c;
  
  Mat x = A.inv() * C;
  Point2i P;
  double ty = x.at<double>(0,0) + 0.5;
  P.y =(int) floor(ty);
  double tx = x.at<double>(1,0) + 0.5;
  P.x = (int) floor(tx);
  //printf("Point P x=%d y=%d\n",P.x,P.y);
  return(P);
}


double partCoveredByTwoParallelLine(Vec4i line1, Vec4i line2)
{
  double frac;
  if(LineType(line1) != LineType(line2))
  {
    printf("ERROR: In Function partCoveredByTwoParallelLine ==> two input lines are not parallel\n");
    exit(0);
  }
  int type;
  type = LineType(line1);
  double angle = FindAngleofaLine(line1);
  angle = abs(angle);
  Point2i P1,Q1,P2,Q2,P12,Q12,P21,Q21;
  
  double CommomPart;
  double UnionPart;

  
  if(type == 1) // Horizontal Line
  {
    P1.x = min(line1[0],line1[2]); 
    if(P1.x == line1[0])
    {
      P1.y = line1[1];
      Q1.x = line1[2];
      Q1.y = line1[3];
    }
    else
    {
      P1.y = line1[3];
      Q1.x = line1[0];
      Q1.y = line1[1];
    }
    
    P2.x = min(line2[0],line2[2]); 
    if(P2.x == line2[0])
    {
      P2.y = line2[1];
      Q2.x = line2[2];
      Q2.y = line2[3];
    }
    else
    {
      P2.y = line2[3];
      Q2.x = line2[0];
      Q2.y = line2[1];
    }
    
    if(angle < 0.5) // Parallel to x axis
    {
      // calculating common part
      int x1 = max(P1.x,P2.x);
      int x2 = min(Q1.x,Q2.x);
      CommomPart = abs((x2-x1)*1.0);
      
      //calculating union part
      x1 = min(P1.x,P2.x);
      x2 = max(Q1.x,Q2.x);
      UnionPart = abs((x2 - x1)*1.0);
      
      if(CommomPart <= 0.0) // No Common Part
      {
	CommomPart = 0.0;
      }
      
      frac = CommomPart/UnionPart;
    //  printf("Part Of Horizontal Line %lf\n",frac);
      return(frac);
    }
    else // If Not parallel to x axis
    {
      linesi temp1 = GetLinefrom2Point(line1);
      linesi temp2 = GetLinefrom2Point(line2);
      linesi P1_perp = GetPerpendicularLinePassingThroughaPoint(temp1,P1);
      linesi Q1_perp = GetPerpendicularLinePassingThroughaPoint(temp1,Q1);
      linesi P2_perp = GetPerpendicularLinePassingThroughaPoint(temp2,P2);
      linesi Q2_perp = GetPerpendicularLinePassingThroughaPoint(temp2,Q2);
      P12 = GetIntersectingPoint(P1_perp,temp2);
      Q12 = GetIntersectingPoint(Q1_perp,temp2);
 //     P21 = GetIntersectingPoint(P2_perp,temp1);
 //     Q21 = GetIntersectingPoint(Q2_perp,temp1);
   //   printf("Line2=>P(%d,%d)Q(%d,%d) Line1=P(%d,%d)Q(%d,%d)\n",line2[0],line2[1],line2[2],line2[3],line1[0],line1[1],line1[2],line1[3]);
    //  printf("Line2=>P(%d,%d)Q(%d,%d) Line12=P(%d,%d)Q(%d,%d)\n",P2.x,P2.y,Q2.x,Q2.y,P12.x,P12.y,Q12.x,Q12.y);
      Point2i  x1,x2;
      x1.x = max(P2.x,P12.x);
      if(x1.x == P2.x )
	x1.y = P2.y;
      else
	x1.y = P12.y;
      
      x2.x = min(Q2.x,Q12.x);
      if(x2.x == Q2.x )
	x2.y = Q2.y;
      else
	x2.y = Q12.y;
      
      CommomPart = DistanceBeteenTwoPoint(x1,x2);
      
      
      x1.x = min(P2.x,P12.x);
      if(x1.x == P2.x )
	x1.y = P2.y;
      else
	x1.y = P12.y;
      
      x2.x = max(Q2.x,Q12.x);
      if(x2.x == Q2.x )
	x2.y = Q2.y;
      else
	x2.y = Q12.y;
      
      UnionPart = DistanceBeteenTwoPoint(x1,x2);
      double d1,d2;
      d1 = DistanceBeteenTwoPoint(P12,Q12);
      d2 = DistanceBeteenTwoPoint(P2,Q2);
      
      
      
      if(UnionPart > (d1+d2))
      {
	CommomPart = 0.0;
      }
      frac = CommomPart/UnionPart;
      
    //  printf("CommomPart = %lf UnionPart = %lf and Part = %lf\n",CommomPart,UnionPart,frac);
      
      return(frac);
    }

  }
  
  
  if(type == 2) // Vertical line
  {
    P1.y = min(line1[1],line1[3]); 
    if(P1.y == line1[1])
    {
      P1.x = line1[0];
      Q1.x = line1[2];
      Q1.y = line1[3];
    }
    else
    {
      P1.x = line1[2];
      Q1.x = line1[0];
      Q1.y = line1[1];
    }
    
    P2.y = min(line2[1],line2[3]); 
    if(P2.y == line2[1])
    {
      P2.x = line2[0];
      Q2.x = line2[2];
      Q2.y = line2[3];
    }
    else
    {
      P2.x = line2[2];
      Q2.x = line2[0];
      Q2.y = line2[1];
    }

    
    if(angle > 89.5) // parallel to y axis
    {
      // calculating common part
      int y1 = max(P1.y,P2.y);
      int y2 = min(Q1.y,Q2.y);
      CommomPart = abs((y2-y1)*1.0);
      
      //calculating union part
      y1 = min(P1.y,P2.y);
      y2 = max(Q1.y,Q2.y);
      UnionPart = abs((y2 - y1)*1.0);
      
      if(CommomPart <= 0.0) // No Common Part
      {
	CommomPart = 0.0;
      }
      
      frac = CommomPart/UnionPart;
    //  printf("Part Of Vertical Line %lf\n",frac);
      return(frac);
    }
    else
    {
      linesi temp1 = GetLinefrom2Point(line1);
      linesi temp2 = GetLinefrom2Point(line2);
      linesi P1_perp = GetPerpendicularLinePassingThroughaPoint(temp1,P1);
      linesi Q1_perp = GetPerpendicularLinePassingThroughaPoint(temp1,Q1);
      linesi P2_perp = GetPerpendicularLinePassingThroughaPoint(temp2,P2);
      linesi Q2_perp = GetPerpendicularLinePassingThroughaPoint(temp2,Q2);
      P12 = GetIntersectingPoint(P1_perp,temp2);
      Q12 = GetIntersectingPoint(Q1_perp,temp2);
      P21 = GetIntersectingPoint(P2_perp,temp1);
      Q21 = GetIntersectingPoint(Q2_perp,temp1);
      
      
     // printf("Line2=>P(%d,%d)Q(%d,%d) Line12=P(%d,%d)Q(%d,%d)\n",P2.x,P2.y,Q2.x,Q2.y,P12.x,P12.y,Q12.x,Q12.y);
      
      Point2i  y1,y2;
      y1.y = max(P2.y,P12.y);
      if(y1.y == P2.y )
	y1.x = P2.x;
      else
	y1.x = P12.x;
      
      y2.y = min(Q2.y,Q12.y);
      if(y2.y == Q2.y )
	y2.x = Q2.x;
      else
	y2.x = Q12.x;
      
      CommomPart = DistanceBeteenTwoPoint(y1,y2);
      
      y1.y = min(P2.y,P12.y);
      if(y1.y == P2.y )
	y1.x = P2.x;
      else
	y1.x = P12.x;
      
      y2.y = max(Q2.y,Q12.y);
      if(y2.y == Q2.y )
	y2.x = Q2.x;
      else
	y2.x = Q12.x;
      
      UnionPart = DistanceBeteenTwoPoint(y1,y2);
      double d1,d2;
      d1 = DistanceBeteenTwoPoint(P12,Q12);
      d2 = DistanceBeteenTwoPoint(P2,Q2);
      
      if(UnionPart > (d1+d2))
      {
	CommomPart = 0.0;
      }
      frac = CommomPart/UnionPart;
      
      
     //  printf("CommomPart = %lf UnionPart = %lf and Part = %lf\n",CommomPart,UnionPart,frac);
      
      return(frac);
    }
  }

}

/*
 * return 1 if inside
 * return 0 if outside
 */
int HLineRectangleTest(Rect R, Vec4i Line, int th)
{
  Point2i left,right;
  left.x = min(Line[0],Line[2]); 
  left.y = max(Line[1],Line[3]);
  right.x = max(Line[0],Line[2]);
  right.y = max(Line[1],Line[3]);
  
  if(PointRectangleTest(R,left)==1 && PointRectangleTest(R,right)==1)
  {
    float ratio = (right.x-left.x)*1.0;
    ratio = (R.width*1.0)/ratio;
    if(ratio > 0.6)
      return 1;
    else
      return 0;
  }
  else if(PointRectangleTest(R,left)==0 && PointRectangleTest(R,right)==0)
  {
    if(left.y >= R.y && left.y <= R.y+R.height)
    {
      float ratio = (right.x-left.x)*1.0;
      ratio = (R.width*1.0)/ratio;
      if(ratio > 0.6)
	return 1;
      else
	return 0;
    }
    else
    {
      Vec4i Line1; 
      Line1[0] = R.x; Line1[1] = R.y; 
      Line1[2] = R.x+R.width; Line1[3] = R.y;
      if(DistanceBetweentwoParallelLine(Line,Line1) <= 2*th)
	return 1;
      else
	return 0;
    }
  }
  else
  {
    if(PointRectangleTest(R,left)==1) // left inside right outside
    {
      float part = (R.x+R.width - left.x)*1.0;
      float maxL = max(R.width*1.0,(right.x-left.x)*1.0);
      part = part/maxL;
      if(part > 0.6)
	return 1;
      else
	return 0;
    }
    else // right is inside left outside
    {
      float part = (right.x - R.x)*1.0;
      float maxL = max(R.width*1.0,(right.x-left.x)*1.0);
      part = part/maxL;
      if(part > 0.6)
	return 1;
      else
	return 0;
    }
  }

}

Mat GetROIImage(Mat image, Rect ROI)
{
  Mat NewImage = Mat(ROI.height+1,ROI.width+1,image.type());
  int r,c;
  r = 0;
  for(int i=ROI.y;i<=ROI.y+ROI.height;i++)
  {
    c = 0;
    for(int j=ROI.x;j<=ROI.x+ROI.width;j++)
    {
      NewImage.at<uchar>(r,c) = image.at<uchar>(i,j);
      c++;
    }
    r++;
  }
 // cout << "OLD Image row and col are " << image.rows << ", " << image.cols << endl;
 // cout << "New Image row and col are " << NewImage.rows << ", " << NewImage.cols << endl;
  return NewImage;
}

vector<Rect> FindBoundRects(Mat BinaryImage, Rect R)
{
  Mat TempImage = FindImageInverse(BinaryImage);
  
  vector<vector<Point> > contours;
  vector<Vec4i> hierarchy;
  
  /// Find contours
  findContours( TempImage, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE, Point(0, 0) );
  TempImage.release();
  
  /// Approximate contours to polygons + get bounding rects and circles
  vector<Rect> boundRect;
  
  for( int j = 0; j < contours.size(); j++ )
  {      
    if(hierarchy[j][3] == -1)
    {
      Rect TR = boundingRect( Mat(contours[j]) );
      if(FindOverlappingRectangles(R,TR)!=0)
	boundRect.push_back(TR);
    }
  }
  
  return boundRect;
}



vector<int> TableHProjProfileTop(Mat image, int y, int x_start, int x_end, int th)
{
 // cout << "IN H Proj up " << endl;
  Rect ROI(0,max(30,5*th),image.cols-1,y-max(30,5*th));
  //Mat fRec = GetROIImage(image,ROI);
  //namedWindow("roi",CV_WINDOW_KEEPRATIO);
  //imshow("roi",fRec);
  //waitKey(100);
  //vector<Rect> boundRect = FindBoundRects(fRec);
  vector<Rect> boundRect = FindBoundRects(image,ROI);
  vector<int> HPP;
  int zerosum_cnt = 0;
  int allsum_cnt = 0;
  int cnt = 0;
  while(y>=0+30)
  {
    int sum = 0;   
    for(int i=x_start;i<=x_end;)
    {
      if(i+th <= x_end)
      {
	if(image.at<uchar>(y,i)==0 && image.at<uchar>(y,i+th)==255)
	{
	  int r;
	  for(r=0;r<boundRect.size();r++)
	  {
	    if(PointRectangleTest(boundRect[r],Point(i,y))==1)
	    {
	      if(boundRect[r].x > (x_start-th) && (boundRect[r].x+ boundRect[r].width) < (x_end+th))
	      {
		sum++; i++;		
	      }
	      else
		i++;
	      break;
	    }
	  }
	  if(r==boundRect.size())
	    i++;
	}
	else if(image.at<uchar>(y,i)==0 && image.at<uchar>(y,i+th)==0)
	{
	  int r1,r2;
	  int r;
	  for(r=0;r<boundRect.size();r++)
	  {
	    if(PointRectangleTest(boundRect[r],Point(i,y))==1)
	    {
	      break;
	    }
	  }
	  r1= r;
	  for(r=0;r<boundRect.size();r++)
	  {
	    if(PointRectangleTest(boundRect[r],Point(i+th,y))==1)
	    {
	      break;
	    }
	  }
	  r2 = r;
	  if(r1< boundRect.size() && r2 < boundRect.size())
	  {
	    if(boundRect[r1].x > (x_start-th) && (boundRect[r1].x+ boundRect[r1].width) < (x_end+th))
	    {
	      if(boundRect[r2].x > (x_start-th) && (boundRect[r2].x+ boundRect[r2].width) < (x_end+th))
	      {
		sum = sum + th;
		i = i + th;
	      }
	      else
	      {
		sum++; i++;
	      }
	    }
	    else
	      i++;
	  }
	  else
	    i++;
	}
	else
	  i++;
      }
      else
      {
	if(image.at<uchar>(y,i)==0)
	{
	  int r;
	  for(r=0;r<boundRect.size();r++)
	  {
	    if(PointRectangleTest(boundRect[r],Point(i,y))==1)
	    {
	      if(boundRect[r].x > (x_start-th) && (boundRect[r].x+ boundRect[r].width) < (x_end+th))
	      {
		sum++;		
	      }
	      break;
	    }
	  }
	}
	i++;
      }
    }
    HPP.push_back(sum);
    //cout << "sum=" << sum << endl;
     if( sum <= (x_end-x_start)*0.05)
    {
     // cout << "here zero_th= "<< 5*th << "zerosum= " << zerosum_cnt  << "and HPP size=" << HPP.size() << endl;
      zerosum_cnt++;
      allsum_cnt = 0;
    }
    else if(sum >= (x_end-x_start)*0.85)
    {
      allsum_cnt++;
      zerosum_cnt = 0;
    }
    else
    {
      zerosum_cnt = 0;
      allsum_cnt = 0;
    }
    if(zerosum_cnt > 2.5*th || allsum_cnt > 4*th)
    {
      if(zerosum_cnt>0)
      {
	for(int m=0;m<=zerosum_cnt;m++)
	{
	  if(!HPP.empty())
	    HPP.pop_back();
	  else
	    break;
	}
      }
      else
      {
	for(int m=0;m<=allsum_cnt;m++)
	{
	  if(!HPP.empty())
	    HPP.pop_back();
	  else
	    break;
	}
      }
      break;
    }
    if(cnt > 5 * th)
    {
      break;
    }
    y--;
    cnt++;
  }
  
  if(HPP.size()==zerosum_cnt || HPP.size()==allsum_cnt)
    HPP.clear();
/*  if(!HPP.empty())
  {
  //  cout << "Resizing happening in Upper part" << endl;
    if(zerosum_cnt>0)
      HPP.resize(HPP.size()-zerosum_cnt);
    else
      HPP.resize(HPP.size()-allsum_cnt);
  }*/
  
  for(int i=0;i<HPP.size();i++)
  {
    if(HPP[i] < 2*th)
      HPP[i] = 0;
  }
  
  return HPP;
  
}



vector<int> TableHProjProfileLow(Mat image, int y, int x_start, int x_end, int th)
{
 // cout << "IN H Proj Low " << endl;
  //cout << "y=" << y << " height= " << image.rows-y-1 << endl;
  Rect ROI(0,y,image.cols-1,image.rows-y-1);
  //Mat fRec = GetROIImage(image,ROI);
  //namedWindow("roi",CV_WINDOW_KEEPRATIO);
  //imshow("roi",fRec);
  //waitKey(100);
  //vector<Rect> boundRect = FindBoundRects(fRec);
  vector<Rect> boundRect = FindBoundRects(image,ROI);
  /*Mat draw;
  //image.copyTo(draw);
  cvtColor(image,draw,CV_GRAY2BGR);
  for(int r=0;r<boundRect.size();r++)
  {
    rectangle(draw,boundRect[r],Scalar(255,0,0),2);
  }
  namedWindow( "blocks", CV_WINDOW_KEEPRATIO );
  imshow( "blocks", draw );
  waitKey(0);*/
  vector<int> HPP;
  int zerosum_cnt = 0;
  int allsum_cnt = 0;
  while(y<image.rows-30)
  {
   // cout << "y=" << y;
    int sum = 0;
    
    for(int i=x_start;i<=x_end;)
    {
      if(i+th <= x_end)
      {
	//cout << "If i+th<x_end" << endl;
	if(image.at<uchar>(y,i)==0 && image.at<uchar>(y,i+th)==255)
	{
	  //cout << "In 1st if" << endl;
	  int r = 0;
	  for(r=0;r<boundRect.size();r++)
	  {
	    if(PointRectangleTest(boundRect[r],Point(i,y))==1)
	    {
	      if(boundRect[r].x > (x_start-th) && (boundRect[r].x+ boundRect[r].width) < (x_end+th))
	      {
		sum++; i++;		
	      }
	      else
		i++;
	      break;
	    }
	  }
	  if(r==boundRect.size())
	    i++;
	}
	else if(image.at<uchar>(y,i)==0 && image.at<uchar>(y,i+th)==0)
	{
	   // cout << "In 2nd if" << endl;
	  int r1,r2;
	  int r;
	  for(r=0;r<boundRect.size();r++)
	  {
	    if(PointRectangleTest(boundRect[r],Point(i,y))==1)
	    {
	      break;
	    }
	  }
	  r1 = r;
	  for(int r=0;r<boundRect.size();r++)
	  {
	    if(PointRectangleTest(boundRect[r],Point(i+th,y))==1)
	    {
	      break;
	    }
	  }
	  r2 = r;
	  if(r1 < boundRect.size() && r2 < boundRect.size())
	  {  
	    if(boundRect[r1].x > (x_start-th) && (boundRect[r1].x+ boundRect[r1].width) < (x_end+th))
	    {
	      if(boundRect[r2].x > (x_start-th) && (boundRect[r2].x+ boundRect[r2].width) < (x_end+th))
	      {
		sum = sum + th;
		i = i + th;
	      }
	      else
	      {
		sum++; i++;
	      }
	    }
	    else
	      i++;
	  }
	  else
	    i++;
	}
	else
	  i++;
      }
      else
      {
	//cout << "else part of if i+th<x_end" << endl;
	if(image.at<uchar>(y,i)==0)
	  sum++;
	i++;
      }
      //cout << " x=" << i;
    }
   // cout << endl;
    //cout << "sum=" << sum  << "for y = " << y << endl;
    HPP.push_back(sum);
   // cout << (x_end-x_start)*0.05 << endl;
   // cout << (x_end-x_start)*0.85 << endl;
    if( sum <= (x_end-x_start)*0.05)
    {
      //cout << "here zero_th= "<< 5*th << "zerosum= " << zerosum_cnt  << "and HPP size=" << HPP.size() << endl;
      zerosum_cnt++;
      allsum_cnt = 0;
    }
    else if(sum >= (x_end-x_start)*0.85)
    {
      allsum_cnt++;
      zerosum_cnt = 0;
    }
    else
    {
      zerosum_cnt = 0;
      allsum_cnt = 0;
    }
    if(zerosum_cnt > 5*th || allsum_cnt > 4*th)
    {
     // cout << 5*th << endl;
      if(zerosum_cnt > 5*th)
      {
	for(int m=0;m<=zerosum_cnt;m++)
	{
	  if(!HPP.empty())
	    HPP.pop_back();
	  else
	    break;
	}
      }
      if(allsum_cnt > 4*th)
      {
	for(int m=0;m<=allsum_cnt;m++)
	{
	  if(!HPP.empty())
	    HPP.pop_back();
	  else
	    break;
	}
      }
      break;
    }
   // if(allsum_cnt > th)
     // break;
    y++;
  }
  
  if(HPP.size()==zerosum_cnt || HPP.size()==allsum_cnt)
    HPP.clear();
  /*if(!HPP.empty())
  {
   // cout << "Resizing happening in bottom" << endl;
    if(zerosum_cnt>0)
      HPP.resize(HPP.size()-zerosum_cnt);
    else
      HPP.resize(HPP.size()-allsum_cnt);
  }*/
  
  for(int i=0;i<HPP.size();i++)
  {
    if(HPP[i] < 2*th)
      HPP[i] = 0;
  }
  
  return HPP;
}


vector<int> TableVerticalProjection(Mat image, int y_start, int y_end, int x_start, int x_end, int th)
{
  //cout << "IN V Proj " << endl;
  Rect ROI(0,y_start,image.cols-1,y_end-y_start);
  //Mat fRec = GetROIImage(image,ROI);
  //namedWindow("roi",CV_WINDOW_KEEPRATIO);
  //imshow("roi",fRec);
  //waitKey(100);
  //vector<Rect> boundRect = FindBoundRects(fRec);
  vector<Rect> boundRect = FindBoundRects(image,ROI);
  //vector<int> VG(x_end-x_start,0);
 // vector<int> VG;
  vector<int> Gaps;
  int zerosum_cnt = 0;
  int nonzerosum_cnt = 0;
  for(int x=x_start;x<=x_end;x++)
  {
    int cnt  = 0; 
    for(int y=y_start;y<=y_end;y++)
    {
      if(image.at<uchar>(y,x)==0)
      {
	for(int r=0;r<boundRect.size();r++)
	  {
	    if(PointRectangleTest(boundRect[r],Point(x,y))==1)
	    {
	      if(boundRect[r].x > (x_start-th) && (boundRect[r].x+ boundRect[r].width) < (x_end+th))
	      {
		cnt++;		
	      }
	      break;
	    }
	  }
      }
	//VG[x_posi]++;
    }
    int val1 = (y_end-y_start)/(2.5*th);
    val1 = 0.5*th*val1;
    int val2 = 2.5*th;
   // cout << "In Gap threshold: Two values are " << val1 << " and " << val2 << endl;
    int zeroth = min(val1,val2);
    if(cnt < zeroth)
    {
     // VG.push_back(0);
      cnt = 0;
    }
    //else
     // VG.push_back(cnt);
    if(cnt == 0)
    {
      nonzerosum_cnt = 0;
      zerosum_cnt++;
      if(zerosum_cnt == 1)
	Gaps.push_back(1);
      else
	Gaps[Gaps.size()-1]++;      
    }
    else
    {
      zerosum_cnt = 0;     
      nonzerosum_cnt--;
      if(nonzerosum_cnt == -1)
	Gaps.push_back(-1);
      else
	Gaps[Gaps.size()-1]--;
    }
  }
  
  cout << "size of Gaps is " << Gaps.size() << endl;
  cout << "Lengths of Gaps are:" << endl;
  
  for(int i=0;i<Gaps.size();i++)
    cout << Gaps[i] << "\t" ;
  cout << endl;
  
  for(int i=0;i<Gaps.size();)
  {
    if(Gaps[i]<0)
    {
      //cout << Gaps[i] << " and th=" << (th)*-1 << endl;
      if(Gaps[i] > (th/2)*-1)
      {
	int temp = Gaps[i];
	Gaps.erase(Gaps.begin()+i);
	if(i-1>=0 && i<Gaps.size())
	{
	  //cout << i << endl;
	  Gaps[i-1] = Gaps[i-1] + Gaps[i] + temp;
	  Gaps.erase(Gaps.begin()+i);
	}
      }
      else
	i++;
    }
    else
    {
      if(Gaps[i] < th/2)
      {
	int temp = Gaps[i];
	Gaps.erase(Gaps.begin()+i);
	if(i-1>=0 && i<Gaps.size())
	{
	  //cout << i << endl;
	  Gaps[i-1] = Gaps[i-1] + Gaps[i] + temp;
	  Gaps.erase(Gaps.begin()+i);
	}
      }
      else
	i++;
    }
  }
  
  for(int i=0;i<Gaps.size();i++)
    cout << Gaps[i] << "\t" ;
  cout << endl;
  //exit(0);
  return Gaps;
  //return VG;
}


Rect GetInterSectionOf2Rect(Rect R1, Rect R2)
{
    Point2i TL1 = R1.tl(); Point2i TL2 = R2.tl();
    Point2i BR1 = R1.br(); Point2i BR2 = R2.br();
    Point2i TL,BR;
    TL.x = max(TL1.x,TL2.x); TL.y = max(TL1.y,TL2.y);
    BR.x = min(BR1.x,BR2.x); BR.y = min(BR1.y,BR2.y);
    Rect IR(TL,BR);
    return IR;
}

Rect GetUnionOf2Rect(Rect R1, Rect R2)
{
    Point2i TL1 = R1.tl(); Point2i TL2 = R2.tl();
    Point2i BR1 = R1.br(); Point2i BR2 = R2.br();
    Point2i TL,BR;
    TL.x = min(TL1.x,TL2.x); TL.y = min(TL1.y,TL2.y);
    BR.x = max(BR1.x,BR2.x); BR.y = max(BR1.y,BR2.y);
    Rect UR(TL,BR);
    return UR;
}

int GetHDistBetween2Rect(Rect R1, Rect R2)
{
    int dist;
    Point2i TL1 = R1.tl(); Point2i TL2 = R2.tl();
    Point2i BR1 = R1.br(); Point2i BR2 = R2.br();
    int val = FindOverlappingRectangles(R1,R2);
    if(val == 0)
    {
      bool inside = false; 
      if(TL1.x >= TL2.x && TL1.x <= BR2.x)
	inside = true;
      if(BR1.x >= TL2.x && BR1.x <= BR2.x)
	inside = true;
      if(inside)
      {
	if(TL1.y > TL2.y)
	  dist = TL1.y - BR2.y;
	else
	  dist = TL2.y - BR1.y;
      }
      else 
	dist = 99999;
    }
    else
      dist = 0;
    
    return dist;
}

typedef struct TableRowsCols
{
  int s;
  int e;
  int zp_cnt;
}TabRC;


float DistancebetweenTwoRowofTable(vector<vector<float> > r1, vector<vector<float> > r2)
{
  //cout << "In row dist " << r1[0].size() << endl;
  vector<float> s1(r1[0].size(),0.0);
  vector<float> s2(r1[0].size(),0.0);
  for(int i=0;i<r1.size();i++)
  {
    for(int j=0;j<r1[i].size();j++)
    {
      s1[j] = s1[j] + r1[i][j];
    }
  }
  
  for(int i=0;i<r2.size();i++)
  {
    for(int j=0;j<r2[i].size();j++)
    {
      s2[j] = s2[j] + r2[i][j];
    }
  }
  
  float dist = 0.0;
  for(int i=0;i<s1.size();i++)
  {
    dist = dist + (s1[i]-s2[i])*(s1[i]-s2[i]);
  }
  dist = sqrt(dist);
  
  return dist;
}

bool GetStructureofaTable(Mat image, Rect &R,int th)
{
  
   Mat SmoothedImg;
  image.copyTo(SmoothedImg);
  SmoothedImg = horizontal_gapfilling(SmoothedImg,th*3);
  Mat tempimg;
  SmoothedImg.copyTo(tempimg);
  cvtColor(tempimg,tempimg,CV_GRAY2BGR);
  
  vector<Rect> NR = FindBoundRects(SmoothedImg,R);
  vector<vector<Rect> > NRowWithCol;
  cout << "NR size is " << NR.size() << endl;
  vector<bool> NR_flag(NR.size(),true);
  vector<Rect> NewRows;
  Mat NewImg = Mat(image.rows,image.cols,image.type(),Scalar(255));
  for(int i=0;i<NR.size();i++)
  {
    //rectangle(tempimg,NR[i],Scalar(0,255,0));
    if(NR_flag[i])
    {
      int rxs,rxe;
      rxs = max(NR[i].x,R.x);
      rxe = min(NR[i].x+NR[i].width,R.x+R.width);
      float partx = ((rxe-rxs)*1.0)/NR[i].width;
      if(partx<0.8)
      {
	NR_flag[i] = false;
      }
      else
      {
	R.x = min(NR[i].x,R.x);
	R.width = max(NR[i].x+NR[i].width,R.x+R.width) - min(NR[i].x,R.x);
	//rectangle(tempimg,R,Scalar(0,0,255));
	for(int r=NR[i].y;r<=NR[i].height;r++)
	{
	  for(int c=NR[i].x;c<=NR[i].width;c++)
	    NewImg.at<uchar>(r,c) = image.at<uchar>(r,c);
	}
      }
    }
  }
  NewImg.copySize(image);
  
  
  
  for(int i=0;i<NR.size();i++)
  {
    if(NR_flag[i])
    {
      //rectangle(tempimg,NR[i],Scalar(0,255,0));
      vector<Rect> ColPerRow;
      Rect R1 = NR[i];
      ColPerRow.push_back(R1);
      for(int j=0;j<NR.size();j++)
      {
	if(i!=j && NR_flag[j])
	{
	  Rect R2 = NR[j];
	  //cout<< "For Rect " << i << " and " << j << endl;
	  //cout << "R1.y=" << R1.y << " R1.y+R1.height=" << R1.y+R1.height << "\t";
	  //cout << "R2.y=" << R2.y << " R2.y+R2.height=" << R2.y+R2.height << endl;
	  if((R1.y>=R2.y && R1.y<=R2.y+R2.height)||(R2.y>=R1.y && R2.y<=R1.y+R1.height)) // same row
	  {
	    /*int c_y = min(RR.y+RR.height,R2.y+R2.height) - max(RR.y,R2.y);
	    if(c_y < th/3)
	    {
	      ColPerRow.push_back(R2);
	      RR = GetUnionOf2Rect(RR,R2);
	      NR_flag[j] =false;
	    }
	    */
	   // cout<< "I am in for Rect " << i << " " << j << endl;
	    Rect RR = GetUnionOf2Rect(R1,R2);
	    NR_flag[j] =false;
	    NR[i] = RR;
	    R1 = RR;
	  }
	}
      }
      rectangle(tempimg,R1,Scalar(255,0,0));
      NewRows.push_back(R1);
      NRowWithCol.push_back(ColPerRow);
      ColPerRow.clear();
    }
  }
  //NR_flag.clear();
  
  if(NewRows.empty())
    return false;
  
 /* namedWindow("show",CV_WINDOW_KEEPRATIO);
  imshow("show",tempimg);
  waitKey(0); */
  
  vector<vector<int> > RC_Data(NewRows.size());
  int max_dis_row = 0;
  for(int i=0;i<NewRows.size();i++)
  {
    int c_x;
    RC_Data[i].resize(NewRows.size(),1);
    int cnt_dis_row = 0;
    for(int j=0;j<NewRows.size();j++)
    {
      if(i!=j)
      {
	c_x = min(NewRows[i].x+NewRows[i].width,NewRows[j].x+NewRows[j].width) - max(NewRows[i].x,NewRows[j].x);
	
	if(c_x < R.width*0.05)
	{
	  cout << "0\t";
	  RC_Data[i][j] = 0;
	  cnt_dis_row++;
	}
	else
	  cout << "1\t";
      }
      else
	cout << "0\t";
    }
    cout << endl;
    if(max_dis_row<cnt_dis_row)
      max_dis_row = cnt_dis_row; 
  }
  cout << "Number of rows " << NewRows.size() << endl;
  cout << "max_dis_row =" << max_dis_row << endl;
  float dissimiratity_index = (max_dis_row*1.0)/NewRows.size();
  cout << "Dissimilarity index = " << dissimiratity_index << endl;
  if(dissimiratity_index>=0.6)
    return false;
  
  int small_row_cnt = 0;
  for(int i=0;i<NewRows.size();i++)
  {
    cout << "Row width " << NewRows[i].width << " Table Width " << R.width << " ratio " << (NewRows[i].width*1.0)/R.width << endl;
    if((NewRows[i].width*1.0)/R.width < 0.6)
      small_row_cnt++;
  }
  cout << (small_row_cnt*1.0)/NewRows.size() << endl;
  if((small_row_cnt*1.0)/NewRows.size()>0.6)
    return false;
  

 /*******************************************************************************************************/
 
  NR.clear();
  NRowWithCol.clear();
  RC_Data.clear();
  NR_flag.clear();
  NewRows.clear();
  /*
  NR = FindBoundRects(image,R);
  NR_flag.resize(NR.size(),true);
  int min_row_posi,max_row_posi;
  for(int i=0;i<NR.size();i++)
  {
    if(NR_flag[i])
    {
      vector<Rect> ColPerRow;
      Rect R1 = NR[i];
      //Rect RR = R1;
      ColPerRow.push_back(R1);
      for(int j=0;j<NR.size();j++)
      {
	if(i!=j && NR_flag[j])
	{
	  Rect R2 = NR[j];
	  if((R1.y>=R2.y && R1.y<=R2.y+R2.height)||(R2.y>=R1.y && R2.y<=R1.y+R1.height)) // same row
	  //if(RR.y<=R2.y+R2.height || RR.y+RR.height>=R2.y) // same row
	  {	    
	      ColPerRow.push_back(R2);
	      Rect RR = GetUnionOf2Rect(RR,R2);
	      NR_flag[j] =false;
	      NR[i] = RR;
	      R1 = RR;
	  }
	}
      }
      NewRows.push_back(R1);
      NRowWithCol.push_back(ColPerRow);
      ColPerRow.clear();
    }
  }
  min_row_posi = 0;
  max_row_posi = 0;
  for(int i=0;i<NewRows.size();i++)
  {
    if(NewRows[i].y<NewRows[min_row_posi].y)
      min_row_posi=i;
    if(NewRows[i].y>NewRows[max_row_posi].y)
      max_row_posi = i;
  }
  float part_x_minr = NewRows[min_row_posi].width/R.width;
  float part_x_maxr = NewRows[max_row_posi].width/R.width;
  
  if(part_x_minr<0.6 && NRowWithCol[min_row_posi].size()==1)
  {
    R.y = NewRows[min_row_posi].y + NewRows[min_row_posi].height;
    R.height = R.height - NewRows[min_row_posi].height;
  }
  
  if(part_x_maxr<0.6 && NRowWithCol[max_row_posi].size()==1)
  {
    R.height = R.height - NewRows[max_row_posi].height;
  }*/
  
  /**********************************************************************************************/
  
  
  vector<int> HP(R.height,0);
  int zc = 0;
  int nzc = 0;
  int r = 0;
  vector<TabRC> TR;
  vector<int> RG;
  for(int i=R.y;i<R.y+R.height;i++)
  {
    for(int j=R.x;j<R.x+R.width;j++)
    {
      //cout << i << j << endl;
      if(image.at<uchar>(i,j)==0)
      {
	//cout << "r=" << r << " and HP[r]= " << HP[r] << endl;
	HP[r]++;
      }
    }
    if(HP[r]==0 || HP[r]<th)
    {
      HP[r] = 0;
      if(zc==0)
	zc = 1;
      else
	zc++;
      nzc = 0;
    }
    else
    {
      TabRC Row;
      if(nzc==0)
      {
	//cout << "Number of rows " << TR.size() << endl;
	if(!TR.empty())
	  RG.push_back(zc);
	nzc = 1;
	Row.s = i;
	Row.e = i;
	Row.zp_cnt = nzc;
	TR.push_back(Row);
      }
      else
      {
	TR[TR.size()-1].zp_cnt++;
	TR[TR.size()-1].e = i;
	nzc++;
	//printf("For row %d s=%d e=%d\n",TR.size()-1,TR[TR.size()-1].s,TR[TR.size()-1].e);
      }
      zc = 0;
    }
    r++;
  }
  
  cout << "TR size" << TR.size() << endl;
  
  if(TR.size()<=1)
    return false;
  
  for(int i=0;i<TR.size();)
  {
    if((TR[i].e - TR[i].s) < th/3)
    {
      TR.erase(TR.begin()+i);
    }
    else
      i++;
  }
  /*
  for(int i=0;i<TR.size();i++)
  {
    if((TR[i].e - TR[i].s) > 3.5*th)
      return false;
  }*/
  
  vector<vector<TabRC> > TableCols(TR.size());
  
  zc = 0;
  nzc = 0;
  int max_col = 1;
  for(int i=0;i<TR.size();i++)
  {
    TabRC Row = TR[i];
    vector<int> VP(R.width,0);
    vector<int> CG;
    int c = 0;
    zc = 0;
    nzc = 0;
    vector<TabRC> RowCol;
    for(int m=R.x;m<R.x+R.width;m++)
    {
      for(int n=Row.s;n<=Row.e;n++)
      {
	if(image.at<uchar>(n,m)==0)
	{
	  VP[c]++;
	}
      }
      if(VP[c] == 0 || (TableCols.size()>2 && VP[c]<th/2))
      {
	VP[c] = 0;
	if(zc == 0 )
	  zc = 1;	  
	else
	  zc++;
	nzc = 0;
      }
      else
      {
	TabRC Col;
	if(nzc == 0)
	{
	  if(!TR.empty())
	    CG.push_back(zc);
	  nzc = 1;
	  Col.s = m;
	  Col.zp_cnt = nzc;
	  RowCol.push_back(Col);
	}
	else
	{
	  RowCol[RowCol.size()-1].zp_cnt++;
	  RowCol[RowCol.size()-1].e = m;
	  nzc++;
	}
	zc = 0;
      }
      c++;
    }
    if(max_col<=RowCol.size())
      max_col = RowCol.size();
    TableCols[i] = RowCol;
    //ColGaps[i] = CG;
  }
  
  int single_col_cnt = 0;
  for(int i=0;i<TableCols.size();)
  {
    if(TableCols[i].size() == 1)
      single_col_cnt++;
    
    if(TableCols[i].size() > 0)
      i++;
    else
      TableCols.erase(TableCols.begin()+i);
  }
  
  if(TableCols.size()==single_col_cnt)
    return false;
  
  if(TableCols.empty())
    return false;
 
  
  
 
  /**********col gap and row gap processing******************/
  
  int min_col_posi = TableCols[0][0].s;
  int max_col_posi = TableCols[0][TableCols[0].size()-1].e;
  float mean_col_width = 0.0;
  int col_cnt =0;
  vector<vector<TabRC> > ColGaps(TR.size());
  for(int i=0;i<TableCols.size();i++)
  {
    cout << "In Row " << i << "( " << TR[i].s << "," << TR[i].e << " ) with height " << TR[i].e-TR[i].s << " Number of columns " << TableCols[i].size() << endl;
    int end = R.x;
    for(int j=0; j<TableCols[i].size();j++)
    {
      if(j<TableCols[i].size())
      {
	TabRC Cgap;
	Cgap.s = end;
	Cgap.e = TableCols[i][j].s;
	Cgap.zp_cnt = Cgap.e - Cgap.s; 
	ColGaps[i].push_back(Cgap);
	end = TableCols[i][j].e;
	if(j==TableCols[i].size()-1)
	{
	  Cgap.s = end;
	  Cgap.e = R.x + R.width;
	  Cgap.zp_cnt = Cgap.e - Cgap.s; 
	  ColGaps[i].push_back(Cgap);
	  end = Cgap.e;
	}
	
      }
      //cout << TableCols[i][j].s << "," << TableCols[i][j].e << "\t";
      //cout << "coverage " << (TableCols[i][j].e-TableCols[i][j].s)*1.0/R.width << "\t";
      if(min_col_posi >= TableCols[i][j].s)
	min_col_posi = TableCols[i][j].s;
      if(max_col_posi <= TableCols[i][j].e)
	max_col_posi = TableCols[i][j].e;
      mean_col_width = mean_col_width + ((TableCols[i][j].e - TableCols[i][j].s)*1.0);
      col_cnt++;
    }
    //cout << endl;
  }
  
  mean_col_width = mean_col_width/col_cnt;
  //int no_col_div = (max_col_posi - min_col_posi)/mean_col_width;
 
  
  
  int no_col_div = R.width/th;
  int col_width = th;
  
  vector<vector<vector<float> > > RowColratio(TableCols.size());
  //vector<vector<vector<float> > > RowColGapratio(TableCols.size());
  //vector<vector<vector<float> > > TableRowVec(TableCols.size());
  for(int i=0;i<TableCols.size();i++)
  {
    vector<vector<float> > rc(TableCols[i].size());
    vector<vector<float> > ColGapVec;
    for(int j=0; j<TableCols[i].size();j++)
      {
	vector<float> cc(no_col_div);
	rc[j] =cc;
	ColGapVec.push_back(cc);
      }
      vector<vector<float> > rcg(ColGaps[i].size());
      /*for(int j=0; j<ColGaps[i].size();j++)
      {
	vector<float> cc(no_col_div);
	rcg[j] = cc;
	ColGapVec.push_back(cc);
      }*/
      RowColratio[i] = rc;
      //RowColGapratio[i] = rcg;
      //TableRowVec[i] = ColGapVec;
  }
  vector<TabRC> col_vals(no_col_div);
  for(int k=0;k<no_col_div;k++)
  {
    if(k==0)
      col_vals[k].s = R.x;
    else
      col_vals[k].s = col_vals[k-1].e+1;
    if(k+1==no_col_div)
    col_vals[k].e = R.x+R.width;
    else
      col_vals[k].e = col_vals[k].s + col_width;
    col_vals[k].zp_cnt = col_vals[k].e - col_vals[k].s;
  }
  
  
  int G1_row_no = 0;
  for(int i=0;i<TableCols.size();i++)
  {
    if(TableCols[i].size() > 1)
      G1_row_no++;
   // int cg_posi = 0;
    for(int j=0; j<TableCols[i].size();j++)
    {
      for(int k=0;k<no_col_div;k++)
      {
	if(TableCols[i][j].e < col_vals[k].s || TableCols[i][j].s > col_vals[k].e) // type t1 and t5
	  RowColratio[i][j][k] = 0.0;
	else
	{
	  if(TableCols[i][j].s < col_vals[k].s && TableCols[i][j].e <= col_vals[k].e && TableCols[i][j].e>=col_vals[k].s) // type t2
	    RowColratio[i][j][k] = ((TableCols[i][j].e - col_vals[k].s)*1.0)/col_vals[k].zp_cnt;
	  else if(TableCols[i][j].e > col_vals[k].e && TableCols[i][j].s>= col_vals[k].s && TableCols[i][j].s <= col_vals[k].e) // type t4
	    RowColratio[i][j][k] = ((col_vals[k].e - TableCols[i][j].s)*1.0)/col_vals[k].zp_cnt;
	  else if(TableCols[i][j].s >= col_vals[k].s && TableCols[i][j].e <= col_vals[k].e) // type t3
	    RowColratio[i][j][k] = ((TableCols[i][j].e - TableCols[i][j].s)*1.0)/col_vals[k].zp_cnt;
	  else// type t6 [i][j].s < [k].s && [i][j].e > [k].e
	    RowColratio[i][j][k] = 1.0;
	}
	//TableRowVec[i][cg_posi][k] = RowColratio[i][j][k];
      }
      //cg_posi++;
    }
    // For Gap Part
   /* for(int j=0; j<ColGaps[i].size();j++)
    {
      for(int k=0;k<no_col_div;k++)
      {
	if(ColGaps[i][j].e < col_vals[k].s || ColGaps[i][j].s > col_vals[k].e) // type t1 and t5
	  RowColGapratio[i][j][k] = 0.0;
	else
	{
	  if(ColGaps[i][j].s < col_vals[k].s && ColGaps[i][j].e <= col_vals[k].e && ColGaps[i][j].e>=col_vals[k].s) // type t2
	    RowColGapratio[i][j][k] = ((ColGaps[i][j].e - col_vals[k].s)*1.0)/col_vals[k].zp_cnt;
	  else if(ColGaps[i][j].e > col_vals[k].e && ColGaps[i][j].s>= col_vals[k].s && ColGaps[i][j].s <= col_vals[k].e) // type t4
	    RowColGapratio[i][j][k] = ((col_vals[k].e - ColGaps[i][j].s)*1.0)/col_vals[k].zp_cnt;
	  else if(ColGaps[i][j].s >= col_vals[k].s && ColGaps[i][j].e <= col_vals[k].e) // type t3
	    RowColGapratio[i][j][k] = ((ColGaps[i][j].e - ColGaps[i][j].s)*1.0)/col_vals[k].zp_cnt;
	  else// type t6 [i][j].s < [k].s && [i][j].e > [k].e
	    RowColGapratio[i][j][k] = 1.0;
	}
	TableRowVec[i][cg_posi][k] = RowColGapratio[i][j][k];
      }
      cg_posi++;
    }*/
    
  }
  
  
  
  float max_dist = (no_col_div*1.0);
  max_dist = sqrt(max_dist);
  cout << "Max Possible Distance Between  two rows in this table is " << max_dist << " with col div " << no_col_div << endl;
  float dist_threshold = max_dist*0.07;
  cout << "Distance Threshold is " << dist_threshold << endl;
  if(G1_row_no==2 || R.width < 15*th)
    dist_threshold = max((float)0.51,dist_threshold);
  //dist_threshold = max((float)0.45,dist_threshold);
  cout << "Distance Threshold is " << dist_threshold << endl;
  cout << "distance between each row of the tables are " << endl;
  float similarity_per = 0.0;
  float total_cnt = 0.0;
  vector<vector<int> > Sim_Vec(RowColratio.size());
  int max_sim_row = 0;
  int N_Proper_rows = 0;
  for(int i=0;i< RowColratio.size();i++)
  {
    if(RowColratio[i].size() > 1)
    {
      N_Proper_rows++;
      int sim_cnt = 0;
      for(int j=0;j< RowColratio.size();j++)
      {
	if(RowColratio[j].size() > 1)
	{
	  //cout << "Finding Similarity for " << i << "\t" << j << "is" << endl; 
	  float dist;
	  dist = DistancebetweenTwoRowofTable( RowColratio[i],  RowColratio[j]);
	  
	  dist = dist/max_dist;
	  cout << dist << "\t";
	  if(i!=j)
	  {
	      if(dist<dist_threshold)
		similarity_per = similarity_per + 1.0;	    
	      total_cnt = total_cnt + 1;
	  }
	  if(dist < dist_threshold)
	  {
	    Sim_Vec[i].push_back(j);
	    sim_cnt++;
	  }
	}
      }
      if(max_sim_row<sim_cnt)
	max_sim_row = sim_cnt;
      cout << endl;
      //cout << "distance between each row of t"
    }
  }
  
  //similarity_per = similarity_per/total_cnt;
  
  cout << "Similarity percentage of rows in the Table is " << similarity_per/total_cnt << endl;
  
  cout << "Max no of Rows that are similar are " << max_sim_row << endl;
  
  similarity_per = (max_sim_row*1.0)/N_Proper_rows;
  
  if(RowColratio.size() == 2 && TableCols[1].size() <=1)
    similarity_per = 0.0;
  
  cout << "New Similarity percentage of rows in the Table is " << similarity_per << endl;
  
  bool tflag = true;
  if(similarity_per<0.54)
    tflag = false;
  
  
  
  return tflag;
  
  /*
  
  similarity_per = 0.0;
  for(int i=0;i< RowColGapratio.size();i++)
  {
    int sim_cnt = 0;
    for(int j=0;j< RowColGapratio.size();j++)
    {
	//cout << "Finding Similarity for " << i << "\t" << j << "is" << endl; 
	float dist;
	if(ColGaps[i].size() > 0  && ColGaps[j].size() > 0)
	  dist = DistancebetweenTwoRowofTable( RowColGapratio[i],  RowColGapratio[j]);
	else
	  dist = 99999.0;
	cout << dist << "\t";
	if(i!=j)
	{
	    if(dist<1.5)
	      similarity_per = similarity_per + 1.0;	    
	    total_cnt = total_cnt + 1;
	}
	if(dist < 1.5)
	{
	  Sim_Vec[i].push_back(j);
	  sim_cnt++;
	}
    }
    if(max_sim_row<sim_cnt)
      max_sim_row = sim_cnt;
    cout << endl;
    //cout << "distance between each row of t"
  }
  
  similarity_per = similarity_per/total_cnt;
  
  cout << "Similarity percentage of rows in the Table is " << similarity_per << endl;
  
  similarity_per = (max_sim_row*1.0)/RowColratio.size();
  
  if(RowColratio.size() == 2 && TableCols[1].size() <=1)
    similarity_per = 0.0;
  
  cout << "New Similarity percentage of rows in the Table is " << similarity_per << endl;
  
  //TableRowVec
  
  similarity_per = 0.0;
  for(int i=0;i< TableRowVec.size();i++)
  {
    int sim_cnt = 0;
    for(int j=0;j< TableRowVec.size();j++)
    {
	//cout << "Finding Similarity for " << i << "\t" << j << "is" << endl; 
	float dist;
	if((ColGaps[i].size() > 0  && ColGaps[j].size() > 0) || (TableCols[i].size() > 0  && TableCols[j].size() > 0))
	  dist = DistancebetweenTwoRowofTable( TableRowVec[i],  TableRowVec[j]);
	else
	  dist = 99999.0;
	cout << dist << "\t";
	if(i!=j)
	{
	    if(dist<1.5)
	      similarity_per = similarity_per + 1.0;	    
	    total_cnt = total_cnt + 1;
	}
	if(dist < 1.5)
	{
	  Sim_Vec[i].push_back(j);
	  sim_cnt++;
	}
    }
    if(max_sim_row<sim_cnt)
      max_sim_row = sim_cnt;
    cout << endl;
    //cout << "distance between each row of t"
  }
  
  similarity_per = similarity_per/total_cnt;
  
  cout << "Similarity percentage of rows in the Table is " << similarity_per << endl;
  
  similarity_per = (max_sim_row*1.0)/RowColratio.size();
  
  if(RowColratio.size() == 2 && TableCols[1].size() <=1)
    similarity_per = 0.0;
  
  cout << "New Similarity percentage of rows in the Table is " << similarity_per << endl;
  */
}

int GetVDistBetween2Rect(Rect R1, Rect R2)
{
    int dist;
    Point2i TL1 = R1.tl(); Point2i TL2 = R2.tl();
    Point2i BR1 = R1.br(); Point2i BR2 = R2.br();
    int val = FindOverlappingRectangles(R1,R2);
    if(val == 0)
    {
      bool inside = false; 
      if(TL1.y >= TL2.y && TL1.y <= BR2.y)
	inside = true;
      if(BR1.y >= TL2.y && BR1.y <= BR2.y)
	inside = true;
      if(inside)
      {
	if(TL1.x > TL2.x)
	  dist = TL1.x - BR2.x;
	else
	  dist = TL2.x - BR1.x;
      }
      else
	dist = 99999;
    }
    else
      dist = 0;
    
    return dist;
}


int PixelsInaRect(Mat image, Rect R)
{
  int pix_cnt = 0;
  for(int i=R.y;i<=R.y+R.height;i++)
  {
    for(int j=R.x;j<=R.x+R.width;j++)
    {
      if(image.at<uchar>(i,j)==0)
	pix_cnt++;
    }
  }
  return pix_cnt;
}


vector<Rect> FindBoundRectNotInVecR(Mat image, vector<Rect> VR)
{
  vector<Rect> BoundRects = FindBoundRects(image,Rect(0,0,image.cols-1,image.rows-1));
  for(int i=0;i<VR.size();i++)
  {
    for(int j=0;j<BoundRects.size();)
    {
      if(FindOverlappingRectangles(VR[i],BoundRects[j])!=0)
      {
	BoundRects.erase(BoundRects.begin()+j);
      }
      else
	j++;
    }
  }
  return BoundRects;
}


vector<Rect> FindTablesWithNoHLines(Mat image, vector<Rect> PTableWithLines, int th)
{
  vector<Rect> BoundRects = FindBoundRectNotInVecR(image,PTableWithLines);
  Mat tempimg;
  image.copyTo(tempimg);
  cvtColor(tempimg,tempimg,CV_GRAY2BGR);
  
  for(int i=0;i<BoundRects.size();)
  {
    if(BoundRects[i].width>image.cols/4)
      BoundRects.erase(BoundRects.begin()+i);
    else
    {
      rectangle(tempimg,BoundRects[i],Scalar(255,0,0));
      i++;
    }
  }
  
  /*namedWindow("show",CV_WINDOW_KEEPRATIO);
  imshow("show",tempimg);
  waitKey(0);
  */
  vector<bool> RectFlags(BoundRects.size(),true);
  
   /******************************/
   
   // generating lines
  
  vector<vector<Rect> > Rows;
  vector<Rect> Lines;
  for(int i=0;i<BoundRects.size();i++)
  {
    if(RectFlags[i])
    {
      Rect R1 = BoundRects[i];
      vector<Rect> row;
      row.push_back(R1);
      //row.push_back(BoundRects[i]);
      for(int j=0;j<BoundRects.size();j++)
      {
	if(RectFlags[j] && i!=j)
	{
	  Rect R2 = BoundRects[j];
	  if(((R1.y>=R2.y && R1.y<=R2.y+R2.height)||(R2.y>=R1.y && R2.y<=R1.y+R1.height))) // same row
	  {
	    if(R1.x<R2.x)
	    {
	      int dist = R2.x-(R1.x+R1.width);
	      if(dist < 10*th && dist > 2*th)
	      {
		Rect RR = GetUnionOf2Rect(R1,R2);
		RectFlags[j] =false;
		R1 = RR;
		row.push_back(R2);
	      }
	    }
	    else
	    {
	      int dist = R1.x-(R2.x+R2.width);
	      if(dist < 10*th && dist > 2*th)
	      {
		Rect RR = GetUnionOf2Rect(R1,R2);
		RectFlags[j] =false;
		R1 = RR;
		row.push_back(R2);
	      }
	    }
	    
	  }
	}
      }
      //rectangle(tempimg,R1,Scalar(0,0,255),3);
      if(row.size()>1)
      {
	rectangle(tempimg,R1,Scalar(0,0,255),3);
	Lines.push_back(R1);
	Rows.push_back(row);
	RectFlags[i] = false;
      }
    }
  }
  
  /*namedWindow("show",CV_WINDOW_KEEPRATIO);
  imshow("show",tempimg);
  waitKey(0);
   */
   
  /******************************/
  
  // Reducing lines based on distance threshold horizontally (left allignment)
  
  vector<vector<Rect> > AllignedLines;
  vector<vector<vector<Rect> > > AllignedLinesWithCols;
  vector<bool> TempLineFlag(Lines.size(),true);
  
  for(int i=0;i<Lines.size();i++)
  {
    if(TempLineFlag[i])
    {
      vector<Rect> AL;
      vector<vector<Rect> > ALC;
      AL.push_back(Lines[i]);
      ALC.push_back(Rows[i]);
      Rect R = Lines[i];
      for(int j=0;j<Lines.size();j++)
      {
	if(i!=j && TempLineFlag[j])
	{
	  if(abs(Lines[i].x-Lines[j].x)<2.5*th && abs((Lines[i].x+Lines[i].width)-(Lines[j].x+Lines[j].width))<2.5*th)
	  {
	    R = GetUnionOf2Rect(R,Lines[j]);
	    AL.push_back(Lines[j]);
	    ALC.push_back(Rows[j]);
	    TempLineFlag[j] = false;
	  }
	}
      }
      //rectangle(tempimg,R,Scalar(0,255,0),3);
      AllignedLines.push_back(AL);
      AllignedLinesWithCols.push_back(ALC);
      //Lines[i] = R;
    }
  }
  
  /*namedWindow("show",CV_WINDOW_KEEPRATIO);
  imshow("show",tempimg);
  waitKey(0);*/
  
  Lines.clear();
  Rows.clear();
  TempLineFlag.clear();
  
  /******************************/
  /* 
  for(int k=0;k<AllignedLines.size();k++)
  {
    Lines = AllignedLines[k];
    Rows = AllignedLinesWithCols[k];
    for(int i=0;i<Lines.size();)
    {
      if(Rows[i].size()<2)
      {
	Lines.erase(Lines.begin()+i);
	Rows.erase(Rows.begin()+i);
      }
      else
	i++;
    }
    AllignedLines[k] = Lines;
    AllignedLinesWithCols[k] = Rows;
  }*/
  
  // Removing Allinged single lines
  for(int k=0;k<AllignedLines.size();)
  {
    if(AllignedLines[k].size()<2)
    {
      AllignedLines.erase(AllignedLines.begin()+k);
      AllignedLinesWithCols.erase(AllignedLinesWithCols.begin()+k);
    }
    else
    {    
      k++;
    }
  }
  
 /* namedWindow("show",CV_WINDOW_KEEPRATIO);
  imshow("show",tempimg);
  waitKey(0); */
  
  //sorting lines based on y posi
  for(int k=0;k<AllignedLines.size();k++)
  {
    Lines = AllignedLines[k];
    Rows = AllignedLinesWithCols[k];
    for(int i=0;i<Lines.size();i++)
    {
      int min = i;
      for(int j=i+1;j<Lines.size();j++)
      {
	if(Lines[j].y<Lines[min].y)
	{
	  min = j;
	}
      }
      if(min!=i)
      {
	 Rect TempR = Lines[i];
	 Lines[i] = Lines[min];
	 Lines[min] = TempR;
	 vector<Rect> tempVR = Rows[i];
	 Rows[i] = Rows[min];
	 Rows[min] = tempVR;
      }
    }
    
    //Lines = AllignedLines[k];
      Rect R = Lines[0];
      for(int i=1;i<Lines.size();i++)
      {
	R = GetUnionOf2Rect(R,Lines[i]);
      }
      rectangle(tempimg,R,Scalar(0,255,0),3);
    
    AllignedLines[k] = Lines;
    AllignedLinesWithCols[k] = Rows;
  }
  Lines.clear();
  Rows.clear();
  
  /*namedWindow("show",CV_WINDOW_KEEPRATIO);
  imshow("show",tempimg);
  waitKey(0);*/
  
  /******************************/
  
  // Removing lines based on distance
  
  vector<Rect> NewTables;
  
  for(int k=0;k<AllignedLines.size();k++)
  {
    Lines = AllignedLines[k];
    Rows = AllignedLinesWithCols[k];
    int cnt = 0;
    
    Rect Table = Lines[0];
    for(int i=0;i<Lines.size()-1;)
    {
      if((Lines[i+1].y-(Lines[i].y+Lines[i].height)) < 5*th)
      {
	cout << "i am in for i=" << i << endl;
	Table = GetUnionOf2Rect(Table,Lines[i+1]);
	i = i + 1;
	cnt++;
      }
      else
      {
	if(cnt>2)
	{
	    NewTables.push_back(Table);
	    rectangle(tempimg,Table,Scalar(0,255,255),3);
	}
	cnt  = 0; 
	i = i + 1;
	Table = Lines[i];
      }
    }
    if(cnt > 2)
    {
      rectangle(tempimg,Table,Scalar(0,255,255),3);
      NewTables.push_back(Table);
    }
  }
  
  /*namedWindow("show",CV_WINDOW_KEEPRATIO);
  imshow("show",tempimg);
  waitKey(0);*/
  
  for(int i=0;i<NewTables.size();)
  {
    if(GetStructureofaTable(image,NewTables[i],th))
    {
      i++;
    }
    else
    {
      NewTables.erase(NewTables.begin()+i);
    }
  }
  
  return NewTables;
  
}


bool FineTuneTheTableRegion(Mat image, Rect &Table, vector<Rect> OtherPart_Rect,  vector<vector<Vec4i> > JoinedLine, int th)
{
  
  vector<Vec4i> LinePerRect;

      for(size_t k = 0; k < JoinedLine.size(); k++)
      {
	for(int j=0;j<JoinedLine[k].size();j++)
	{
	  int typei = LineType(JoinedLine[k][j]);
	  if(typei == 1) // Horizontal lines
	  {
	    if(HLineRectangleTest(Table,JoinedLine[k][j], th))
	    {
	      LinePerRect.push_back(JoinedLine[k][j]);
	    }
	  }
	}
      }


      for(int m=0;m<LinePerRect.size();m++)
      {
	int min = m;
	for(int n=m+1;n<LinePerRect.size();n++)
	{
	  if(LinePerRect[n][1]<LinePerRect[min][1])
	  {
	    min = n;
	  }
	}
	if(min!=m)
	{
	  Vec4i L = LinePerRect[m];
	  LinePerRect[m] = LinePerRect[min];
	  LinePerRect[min] = L;
	}
      }
  
  
  
  

      for( int j = 0; j < OtherPart_Rect.size(); j++ )
      {
	if(FindOverlappingRectangles(Table,OtherPart_Rect[j])!=0)
	{
	  if(OtherPart_Rect[j].height>image.rows/6 && OtherPart_Rect[j].width>image.cols/6)
	  {
	    Rect R = Rect(Table.x,OtherPart_Rect[j].y,Table.width,OtherPart_Rect[j].height);
	    Table = R;	    
	  }
	  if(FindOverlappingRectangles(Table,OtherPart_Rect[j])==1)
	  {
	    float ratio = (Table.area()*1.0)/OtherPart_Rect[j].area();
	    if(ratio<0.3)
	      return false;
	  }
	}
      }

  
  vector<Rect> BoundRects = FindBoundRects(image,Table);
  vector<bool> RectFlags(BoundRects.size(),true);
  
   /******************************/
   
   // generating lines
  
  vector<vector<Rect> > Rows;
  vector<Rect> Lines;
  for(int i=0;i<BoundRects.size();i++)
  {
    if(RectFlags[i])
    {
      Rect R1 = BoundRects[i];
      vector<Rect> row;
      row.push_back(R1);
      for(int j=0;j<BoundRects.size();j++)
      {
	if(RectFlags[j] && i!=j)
	{
	  Rect R2 = BoundRects[j];
	  if(((R1.y>=R2.y && R1.y<=R2.y+R2.height)||(R2.y>=R1.y && R2.y<=R1.y+R1.height))) // same row
	  {
	    Rect RR = GetUnionOf2Rect(R1,R2);
	    RectFlags[j] =false;
	    R1 = RR;
	    row.push_back(R2); 
	  }
	}
      }
      Lines.push_back(R1);
      Rows.push_back(row);
      RectFlags[i] = false;
    }
  }
  
  /*namedWindow("show",CV_WINDOW_KEEPRATIO);
  imshow("show",tempimg);
  waitKey(0);
   */
   
  for(int k=0;k<Lines.size();k++)
  {
    vector<Rect> row = Rows[k];
    for(int i=0;i<row.size();i++)
    {
      int min = i;
      for(int j=i+1;j<row.size();j++)
      {
	if(row[j].x<row[min].x)
	{
	  min = j;
	}
      }
      if(min!=i)
      {
	 Rect TempR = row[i];
	 row[i] = row[min];
	 row[min] = TempR;
      }
    }
    Rows[k] = row;
  }
  
  for(int i=0;i<Lines.size();i++)
  {
      int min = i;
      for(int j=i+1;j<Lines.size();j++)
      {
	if(Lines[j].y<Lines[min].y)
	{
	  min = j;
	}
      }
      if(min!=i)
      {
	 Rect TempR = Lines[i];
	 Lines[i] = Lines[min];
	 Lines[min] = TempR;
	 vector<Rect> tempVR = Rows[i];
	 Rows[i] = Rows[min];
	 Rows[min] = tempVR;
      }
  }
  
  // Now we have H st lines in a table in sorted way
  // we have text lines in a table in sorted way
  // we have to check text lines above 1st st line 1) left alligned? 2) posi? 3) length 4) coverage over the table
  // check number of st line in the table
  // if more than 2 
  // check posi of last st line
  // if last st line is near to the end of the table then
  // we have to check text lines above 1st st line 1) left alligned? 2) posi? 3) length 4) coverage over the table
  
  
  //Text Lines above the top most line
  //TLAtl 
  int topi;
  int boti;
  vector<bool> LineFlag(Lines.size(),true);
  if(!LinePerRect.empty())
  {
    cout << "In Top\n";
    if(LinePerRect[0][1] < (Table.y+Table.height*0.2))
    {  
      for(int i=0;i<Lines.size();i++)
      {
	
	if((Lines[i].y+Lines[i].height)<LinePerRect[0][1]) // Text line is above top most horizontal line
	{
	  float lenratio = (Lines[i].width*1.0)/Table.width;
	  float textdensity = 0.0;
	  for(int j=0;j<Rows[i].size();j++)
	  {
	    textdensity = textdensity + Rows[i][j].width;
	  }
	  textdensity = textdensity/Table.width;
	  int strat_point = Lines[i].x;
	  int dist_from_TableLeft = abs(Lines[i].x-Table.x);
	  int end_point = Lines[i].x+Lines[i].width;
	  int dist_from_TableRight = abs((Lines[i].x+Lines[i].width)-(Table.x+Table.width));
	  if((textdensity>0.75 && lenratio>0.75) || (textdensity<0.31 && lenratio<0.31))
	    LineFlag[i] = false;
	  if(dist_from_TableLeft>4*th && dist_from_TableRight>4*th)
	    LineFlag[i] = false;
	  
	  cout << "textdensity=" << textdensity << " lenratio=" << lenratio << " dist_from_TableLeft=" << dist_from_TableLeft << " dist_from_TableRight=" << dist_from_TableRight << endl;
	}
	else
	{
	  topi = i;
	  break;
	}
      }
    }
    if(LinePerRect.size()>1)
    {
      cout << "In Bottom\n";
      // Check position of last line
      if(LinePerRect[LinePerRect.size()-1][1]>(Table.y+Table.height*0.7))
      {
	for(int i=Lines.size()-1;i>0;i--)
	{
	  if(Lines[i].y>LinePerRect[LinePerRect.size()-1][1]) // Text line below last H st line
	  {
	    float lenratio = (Lines[i].width*1.0)/Table.width;
	    float textdensity = 0.0;
	    for(int j=0;j<Rows[i].size();j++)
	    {
	      textdensity = textdensity + Rows[i][j].width;
	    }
	    textdensity = textdensity/Table.width;
	    int strat_point = Lines[i].x;
	    int dist_from_TableLeft = abs(Lines[i].x-Table.x);
	    int end_point = Lines[i].x+Lines[i].width;
	    int dist_from_TableRight = abs((Lines[i].x+Lines[i].width)-(Table.x+Table.width));
	    if((textdensity>0.75 && lenratio>0.65) || (textdensity<0.31 && lenratio<0.31))
	      LineFlag[i] = false;
	    if(dist_from_TableLeft>4*th && dist_from_TableRight>4*th)
	      LineFlag[i] = false;
	    cout << "textdensity=" << textdensity << " lenratio=" << lenratio << " dist_from_TableLeft=" << dist_from_TableLeft << " dist_from_TableRight=" << dist_from_TableRight << endl; 
	  }
	  else
	  {
	    boti = i;
	    break;
	  }
	}
      }
    }
  }
  
  if(topi!=0)
  {
    for(int i=0;i<topi;i++)
    {
      if(!LineFlag[i])
      {
	int y;
	if(i+1!=topi)
	  y = Lines[i+1].y;
	else
	  y = LinePerRect[0][1];
	Table = Rect(Table.x,y,Table.width,Table.height-(y-Table.y));
      }
    }
    
  }
  if(boti!=Lines.size()-1)
  {
    for(int i=Lines.size()-1;i>boti;i--)
    {
      if(!LineFlag[i])
      {
	int y;
	int height;
	if(i-1!=boti)
	{
	  y = Lines[i-1].y+Lines[i-1].height;
	  height = y-Table.y;
	}
	else
	{
	  y = LinePerRect[LinePerRect.size()-1][1];
	  height = y-Table.y;
	}
	Table = Rect(Table.x,Table.y,Table.width,height);
      }
    }
  }
  
  return true;
  
}


void FineTuneTheTableRegion(Mat image, Rect &Table, float th)
{
 // int hgap = (int) ((th*0.5)+0.5);
 // Mat SmoothedImg = horizontal_gapfilling(image,hgap);
 // vector<Rect> BoundRects = FindBoundRects(SmoothedImg,Table);
  vector<Rect> BoundRects = FindBoundRects(image,Table);
  vector<bool> RectFlags(BoundRects.size(),true);
  
   /******************************/
   
   // generating lines
  
  vector<vector<Rect> > Rows;
  vector<Rect> Lines;
  for(int i=0;i<BoundRects.size();i++)
  {
    if(RectFlags[i])
    {
      Rect R1 = BoundRects[i];
      vector<Rect> row;
      row.push_back(R1);
      for(int j=0;j<BoundRects.size();j++)
      {
	if(RectFlags[j] && i!=j)
	{
	  Rect R2 = BoundRects[j];
	  if(((R1.y>=R2.y && R1.y<=R2.y+R2.height)||(R2.y>=R1.y && R2.y<=R1.y+R1.height))) // same row
	  {
	    Rect RR = GetUnionOf2Rect(R1,R2);
	    RectFlags[j] =false;
	    R1 = RR;
	    row.push_back(R2); 
	  }
	}
      }
      Lines.push_back(R1);
      Rows.push_back(row);
      RectFlags[i] = false;
    }
  }
  
  /*namedWindow("show",CV_WINDOW_KEEPRATIO);
  imshow("show",tempimg);
  waitKey(0);
   */
   
  for(int k=0;k<Lines.size();k++)
  {
    vector<Rect> row = Rows[k];
    for(int i=0;i<row.size();i++)
    {
      int min = i;
      for(int j=i+1;j<row.size();j++)
      {
	if(row[j].x<row[min].x)
	{
	  min = j;
	}
      }
      if(min!=i)
      {
	 Rect TempR = row[i];
	 row[i] = row[min];
	 row[min] = TempR;
      }
    }
    Rows[k] = row;
  }
  
  for(int i=0;i<Lines.size();i++)
  {
      int min = i;
      for(int j=i+1;j<Lines.size();j++)
      {
	if(Lines[j].y<Lines[min].y)
	{
	  min = j;
	}
      }
      if(min!=i)
      {
	 Rect TempR = Lines[i];
	 Lines[i] = Lines[min];
	 Lines[min] = TempR;
	 vector<Rect> tempVR = Rows[i];
	 Rows[i] = Rows[min];
	 Rows[min] = tempVR;
      }
  }
  
  
    vector<Rect> row1 = Rows[0];
   // vector<float> gaps1;
    float max_gap = (row1[1].x-(row1[0].x+row1[0].width))*1.0;
    int row1_colcnt = 0;
    for(int i=1;i<row1.size();i++)
    {
      if((row1[i].x-(row1[i-1].x+row1[i-1].width)) > 1.5*th)
	row1_colcnt++;
   //   gaps1.push_back((row1[i].x-(row1[i-1].x+row1[i-1].width))*1.0);
     // if(max_gap<(row1[i].x-(row1[i-1].x+row1[i-1].width))*1.0)
	//max_gap = (row1[i].x-(row1[i-1].x+row1[i-1].width))*1.0;
    }
   // for(int i=0;i<gaps1.size();i++)
     // gaps1[i] = gaps1[i]/max_gap;
    vector<Rect> row2 = Rows[Rows.size()-1];
  //  vector<float> gaps2;
    int row2_colcnt = 0;
    max_gap = (row2[1].x-(row2[0].x+row2[0].width))*1.0;
    for(int i=1;i<row2.size();i++)
    {
      if((row2[i].x-(row2[i-1].x+row2[i-1].width)) > 1.5*th)
	row2_colcnt++;
     // gaps2.push_back((row2[i].x-(row2[i-1].x+row2[i-1].width))*1.0);
    //  if(max_gap<(row2[i].x-(row2[i-1].x+row2[i-1].width))*1.0)
	//max_gap = (row2[i].x-(row2[i-1].x+row2[i-1].width))*1.0;
    }
    
  
  
  /******************************/
  
  // Reducing lines based on distance threshold horizontally (left allignment)
  bool line1_flag = true;
  bool line2_flag = true;
  if(Lines.size()>1)
  {
    if(Lines[0].x > Lines[1].x)
    {
      if(abs(Lines[1].x-Lines[0].x)<4*th || abs((Lines[1].x+Lines[1].width)-(Lines[0].x+Lines[0].width))<4*th)
	line1_flag = true;
      else
      {
	line1_flag = false;
	Table = Rect(Table.x,Lines[1].y,Table.width,Table.height-(Lines[0].height+(Lines[1].y-(Lines[0].height+Lines[0].y))));
      }
    }
    if(Lines[Lines.size()-1].x > Lines[Lines.size()-2].x)
    {
      if(abs(Lines[Lines.size()-1].x-Lines[Lines.size()-2].x)<4*th || abs((Lines[Lines.size()-1].x+Lines[Lines.size()-1].width)-(Lines[Lines.size()-2].x+Lines[Lines.size()-2].width))<4*th)
	line2_flag = true;
      else
      {
	line2_flag = false;
	Table = Rect(Table.x,Table.y,Table.width,Table.height-(Lines[Lines.size()-1].height+(Lines[Lines.size()-1].y-(Lines[Lines.size()-2].height+Lines[Lines.size()-2].y))));
      }
    }
  }
  
  if(line1_flag && row1_colcnt > 0)
  {
    line1_flag = false;
    Table = Rect(Table.x,Lines[1].y,Table.width,Table.height-(Lines[0].height+(Lines[1].y-(Lines[0].height+Lines[0].y))));
  }
  
  if(line2_flag && row2_colcnt > 0)
  {
    line2_flag = false;
    Table = Rect(Table.x,Table.y,Table.width,Table.height-(Lines[Lines.size()-1].height+(Lines[Lines.size()-1].y-(Lines[Lines.size()-2].height+Lines[Lines.size()-2].y))));
  
  }
    
  
  /******************************/
}



int main(int argc, char *argv[])
{
  
  char *substring;
  substring = input_image_name_cut(argv[1]); 
  makedir(substring);
  
  char *name,*output,*tempname;
  
  Mat src = imread(argv[1],1);
  Mat BinaryImage = binarization(src,4);
  
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"BinaryImage.png");
  imwrite(name,BinaryImage);
  
  Mat uniform_background;
  
  uniform_background = foreground_masked_image(src,BinaryImage);
  
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"uniform_background.png");
  imwrite(name,uniform_background);
  Mat ErodedImage;
 /* ErodedImage = dilation(BinaryImage);
  name = CreateNameIntoFolder(substring,"ErodedImageBin.png");
  imwrite(name,ErodedImage);*/
  /*Mat ErodedImage = Erosion(0,1,BinaryImage);
  
  name = CreateNameIntoFolder(substring,"ErodedImage.png");
  imwrite(name,ErodedImage);
  
  Mat TempImage = FindImageInverse(ErodedImage);*/
  Mat TempImage = FindImageInverse(BinaryImage);
  
  vector<vector<Point> > contours;
  vector<Vec4i> hierarchy;
  
  /// Find contours
  findContours( TempImage, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE, Point(0, 0) );
  TempImage.release();
  
  /// Approximate contours to polygons + get bounding rects and circles
  vector<vector<Point> > contours_poly( contours.size() );
  vector<Rect> boundRect( contours.size() );
  vector<RotatedRect> boundRotatedRect( contours.size() );
  
  for( int j = 0; j < contours.size(); j++ )
  {      
    contours_poly[j] = contours[j];
    //approxPolyDP( Mat(contours[j]), contours_poly[j], 3, true );
    boundRect[j] = boundingRect( Mat(contours_poly[j]) );
    boundRotatedRect[j] = minAreaRect( Mat(contours_poly[j]) );
  }
  
  
  cout << (src.cols)/15 << (src.rows)/15 << endl;
  
  src.copyTo(TempImage);
  vector<bool> validBlock(contours.size(),false);
  vector<int> height;
  vector<int> width;
  float avg_height = 0.0;
  for( int j = 0; j< contours.size(); j++ )
  {
    Rect TempRect;
    TempRect = boundRect[j];
    if(!((TempRect.width > (src.cols)/15 || TempRect.height > (src.rows)/15) 
					 ||
	(TempRect.width > min((src.cols)/25,(src.rows)/25) && TempRect.height > min((src.cols)/25,(src.rows)/25))
      ))
    {
      if(hierarchy[j][3] == -1) // 
      {
	validBlock[j] = true;
	height.push_back(TempRect.height);
	width.push_back(TempRect.width);
	rectangle( TempImage, boundRect[j].tl(), boundRect[j].br(), Scalar(255,0,255), 2, 8, 0 );
      }
    }
  }
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"ValidVlock.png");
  
  imwrite(name,TempImage);
  TempImage.release();
  
  double MedianHeight = FindMedian<int>(height);
  double MedianWidth = FindMedian<int>(width);
  double mean_ht = GetVecMean<int>(height);
  cout << "Meadian Height = " << MedianHeight << endl;
  cout << "Mean Height = " << mean_ht << endl;
  
 
  height.clear();
  width.clear();
  
  Mat NewImage = Mat(src.rows,src.cols,src.type(),Scalar(255,255,255));
  Mat NewBinaryDst = Mat(src.rows,src.cols,BinaryImage.type(),Scalar(255));
  
  
  for( int j = 0; j< contours.size(); j++ )
  {
    if(validBlock[j])
    {
      int p = 0 ;
      for(int m=boundRect[j].y;m<boundRect[j].y+boundRect[j].height;m++)
      {
	int q = 0;
	for(int n=boundRect[j].x;n<boundRect[j].x+boundRect[j].width;n++)
	{
	  int temp_col = boundRect[j].width;
	  bool measure_dist;
	  if((pointPolygonTest(contours_poly[j],Point(n,m),measure_dist) >= 0.0) && BinaryImage.data[m*BinaryImage.cols+n]==0)
	  {
	    NewImage.at<Vec3b>(m,n)[0] = uniform_background.at<Vec3b>(m,n)[0];
	    NewImage.at<Vec3b>(m,n)[1] = uniform_background.at<Vec3b>(m,n)[1];
	    NewImage.at<Vec3b>(m,n)[2] = uniform_background.at<Vec3b>(m,n)[2];
	    NewBinaryDst.at<uchar>(m,n) = BinaryImage.at<uchar>(m,n);
	  }
	}
      }
    }
  }
  NewBinaryDst = Erosion(0,1,NewBinaryDst);
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"NewBinaryImage.png");
  imwrite(name, NewBinaryDst);
  
  int hgap,vgap;
  //hgap =(int) ((mean_ht*1.5) + 0.5);
  hgap = (int) ((mean_ht*1.2)+0.5);
  Mat SmoothedImg = horizontal_gapfilling(NewBinaryDst,hgap);
 // SmoothedImg = vertical_gapfilling(SmoothedImg,mean_ht/3);
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"HorizontalBlockImage.png");
  imwrite(name, SmoothedImg);
  /*vgap =(int) ((2.5*mean_ht) + 0.5);
  SmoothedImg = vertical_gapfilling(SmoothedImg,vgap);
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"ProbableTextBlocks.png");
  imwrite(name, SmoothedImg);
  */
  
  TempImage = FindImageInverse(SmoothedImg);
  /// Find contours
  findContours( TempImage, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, Point(0, 0) );	
  
  contours_poly.resize(contours.size());
  boundRect.resize(contours.size());
  boundRotatedRect.resize(contours.size());
  
  vector<int> NewHtData;
  for( int j = 0; j < contours.size(); j++ )
  { 
    approxPolyDP( Mat(contours[j]), contours_poly[j], 3, true );
    boundRect[j] = boundingRect( Mat(contours_poly[j]) );
    boundRotatedRect[j] = minAreaRect( Mat(contours_poly[j]) );
  }
  
  Mat ProbableTableText;
  ProbableTableText = Mat(src.rows,src.cols,CV_8UC1,Scalar(255));
 // src.copyTo(TestImage);
  for( int j = 0; j < contours.size(); j++ )
  {
    if(hierarchy[j][3] == -1)
    {
      Rect R = boundRect[j];
      NewHtData.push_back(R.height);
      if(R.width < (TempImage.cols/3.6) && R.height < mean_ht*3) // removing large lines
      {
	if(R.height >= mean_ht/2) // removing noise
	{
	  //rectangle(TestImage,R.tl(),R.br(),Scalar(255,120,50));
	  int p = 0 ;
	  for(int m=boundRect[j].y;m<boundRect[j].y+boundRect[j].height;m++)
	  {
	    int q = 0;
	    for(int n=boundRect[j].x;n<boundRect[j].x+boundRect[j].width;n++)
	    {
	      int temp_col = boundRect[j].width;
	      bool measure_dist;
	      if((pointPolygonTest(contours_poly[j],Point(n,m),measure_dist) >= 0.0) && BinaryImage.data[m*BinaryImage.cols+n]==0)
	      {	 
		ProbableTableText.at<uchar>(m,n) = BinaryImage.at<uchar>(m,n);
	      }
	    }
	  }
	}
      }
    }
  }
  
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"ProbableTableTextBlocks.png");
  imwrite(name, ProbableTableText);
  
  SmoothedImg = horizontal_gapfilling(ProbableTableText,hgap);
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"ProbableTableTextBlocksHGap.png");
  imwrite(name, SmoothedImg);
  
  Mat OtherPart = Mat(src.rows,src.cols,BinaryImage.type(),Scalar(255));
  
  for(int i=0;i<BinaryImage.rows;i++)
  {
    for(int j=0;j<BinaryImage.cols;j++)
    {
      if(BinaryImage.data[i*BinaryImage.cols+j] == 0 && NewBinaryDst.data[i*BinaryImage.cols+j] != 0)
      {
	OtherPart.data[i*BinaryImage.cols+j] = 0;
      }
    }
  }
  
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"OtherPart.png");
  imwrite(name, OtherPart);
  
  if(!NewHtData.empty())
  {
    mean_ht = GetVecMean<int>(NewHtData);
    mean_ht = max(10.0,(3*mean_ht)/5);
  }
  
  
  cout << "New Height of the Doc is " << mean_ht << endl;
  
  Mat LineImage;
  
 
  /*******************Detect Lines *******************/
  
  
  Mat BoundaryImage = boundaryextraction(OtherPart);
  
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"Boundary_OtherPart.png");
  imwrite(name, BoundaryImage);
  
  
  /*ErodedImage = Erosion(1,2,OtherPart);
  TempImage = FindImageInverse(ErodedImage);*/
  /*ErodedImage = Erosion(1,2,BoundaryImage);
  TempImage = FindImageInverse(ErodedImage);*/
  TempImage = FindImageInverse(BoundaryImage);
 // TempImage = FindImageInverse(OtherPart);
  
  
  cvtColor( TempImage, LineImage, CV_GRAY2BGR );
  vector<Vec4i> lines;
    HoughLinesP( TempImage, lines, 1, CV_PI/180, 80, min(src.cols/25,src.rows/25), 2*mean_ht );
    cout << "Number of lines are " << lines.size() << endl;
    for( size_t i = 0; i < lines.size(); i++ )
    {
        line( LineImage, Point(lines[i][0], lines[i][1]),
            Point(lines[i][2], lines[i][3]), Scalar(0,0,255), 3, 8 );
    }
    TempImage.release();
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"DetectedLine.png");
  imwrite(name, LineImage);
  
  LineImage.release();
  src.copyTo(LineImage);
  
  /*for( size_t i = 0; i< lines.size();  )
  {
    double angle = FindAngleofaLine(lines[i]);
    cout << angle << endl;
    if((angle > -1.0 && angle < .95) || (angle > 88.0 && angle < -88.0))
      i++;
    else
      lines.erase(lines.begin()+i);
  }
  for( size_t i = 0; i< lines.size();  i++)
  {
    double angle = FindAngleofaLine(lines[i]);
    cout << "line%d" << i << "angle=" << angle << endl;
  }*/
  
  vector<vector<Vec4i> > ParallelLines = FindParallelLines(lines,4.0);
  
    for( size_t i = 0; i < ParallelLines.size(); i++ )
    {
      Scalar color = Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) );
      vector<Vec4i> TLines = ParallelLines[i];
      for(size_t j = 0; j < TLines.size(); j++)
      {
        line( LineImage, Point(TLines[j][0], TLines[j][1]),
            Point(TLines[j][2], TLines[j][3]), color, 3, 8 );
      }
    }
    
    name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"ParallelLines.png");
  imwrite(name, LineImage);
  LineImage.release();
  
  
  //vector<vector<vector<Vec4i> > > ConnectedLines = FindConnectedLines(ParallelLines,100);
  //vector<vector<vector<Vec4i> > > ConnectedLines = FindConnectedLines(ParallelLines);
  vector<vector<vector<Vec4i> > > ConnectedLines = FindConnectedLinesNew(ParallelLines,15); // 15 for POD // 30 for ICDAR
  ConnectedLines = SplitConnectedLines(ConnectedLines,30);
  
  src.copyTo(LineImage);
  
  for(size_t k = 0; k < ConnectedLines.size(); k++)
  {
    for( size_t i = 0; i < ConnectedLines[k].size(); i++ )
    {
      Scalar color = Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) );
      vector<Vec4i> TLines = ConnectedLines[k][i];
      for(size_t j = 0; j < TLines.size(); j++)
      {
        line( LineImage, Point(TLines[j][0], TLines[j][1]),
            Point(TLines[j][2], TLines[j][3]), color, 3, 8 );
      }
    }
  }
    
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"ConnectedLines.png");
  imwrite(name, LineImage);
  LineImage.release();
  
  
 // vector<vector<Vec4i> > JoinedLine = JoinedConnectedLine(ConnectedLines);
  vector<vector<Vec4i> > JoinedLine = JoinedConnectedLineNew(ConnectedLines);
  
  
  src.copyTo(LineImage);
  for(size_t k = 0; k < JoinedLine.size(); k++)
  {
    Scalar color = Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) );
    for(int i=0;i<JoinedLine[k].size();i++)
    {
     
      line( LineImage, Point(JoinedLine[k][i][0], JoinedLine[k][i][1]),
	     Point(JoinedLine[k][i][2], JoinedLine[k][i][3]), color, 3, 8 );
    }
  }
  
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"JoinedLine.png");
  imwrite(name, LineImage);
  LineImage.release();
  
  
  /************8 Connected componets of Other Parts ***********************************/
  
  TempImage.release();
  OtherPart.copyTo(TempImage);
  TempImage = horizontal_gapfilling(TempImage,mean_ht/2);
  TempImage = vertical_gapfilling(TempImage,mean_ht/2);
  TempImage = FindImageInverse(TempImage);
  vector<vector<Point> > contours_otherpart;
  vector<Vec4i> hierarchy_otherpart;
  
  
  /// Find contours
  findContours( TempImage, contours_otherpart, hierarchy_otherpart, CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE, Point(0, 0) );
  TempImage.release();
  
  /// Approximate contours to polygons + get bounding rects and circles
  vector<vector<Point> > contours_poly_otherpart( contours_otherpart.size() );
  vector<Rect> boundRect_otherpart( contours_otherpart.size() );
  vector<RotatedRect> boundRotatedRect_otherpart( contours_otherpart.size() );
  
  vector<bool> ValidOtherParts(contours_otherpart.size(),false);
  for( int j = 0; j < contours_otherpart.size(); j++ )
  {      
    contours_poly_otherpart[j] = contours_otherpart[j];
    //approxPolyDP( Mat(contours[j]), contours_poly[j], 3, true );
    boundRect_otherpart[j] = boundingRect( Mat(contours_poly_otherpart[j]) );
   // boundRotatedRect_otherpart[j] = minAreaRect( Mat(contours_poly_otherpart[j]) );
    if(hierarchy_otherpart[j][3] == -1)
      ValidOtherParts[j] = true;
  }
  
  
  for( int j = 0; j < contours_otherpart.size(); j++ )
  {
    if(ValidOtherParts[j])
    {
      cout << "Going in " << j << "block" << endl;
      for( int i = 0; i < contours_otherpart.size(); i++ )
      {
	if(ValidOtherParts[i] && i!=j)
	{
	 // cout << "Going in " << j << "block" << endl;
	  int val = FindOverlappingRectangles(boundRect_otherpart[j],boundRect_otherpart[i]);
	  if(val==1)
	    ValidOtherParts[j] = false;
	  else if(val==2)
	    ValidOtherParts[i] = false;
	}
      }
    }
  }

  src.copyTo(LineImage);
  for( int j = 0; j < contours_otherpart.size(); j++ )
      {
	if(ValidOtherParts[j])
	{
	  rectangle(LineImage,boundRect_otherpart[j],Scalar(120,12,79),3,8,0);
	}
      }
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"bounded_Otherpart.png");
  imwrite(name, LineImage);
  LineImage.release();
  
 /***********************Split joined lines w.r.t. CC of Other Part ******************************/
 
 vector<vector<Vec4i> > NewJoinedLine(JoinedLine.size());
 
 src.copyTo(LineImage);
 
 for(size_t k = 0; k < JoinedLine.size(); k++)
  {
    for(int i=0;i<JoinedLine[k].size();i++)
    {
      Vec4i Line = JoinedLine[k][i];
      bool splitl = false;
      for( int j = 0; j < contours_otherpart.size(); j++ )
      {
	if(ValidOtherParts[j])
	{
	  rectangle(LineImage,boundRect_otherpart[j],Scalar(120,255,79),10,8,0);
	  int P1R1 = NewPointRectangleTest(boundRect_otherpart[j],Point(Line[0],Line[1]),(float)mean_ht);
	  int P2R1 = NewPointRectangleTest(boundRect_otherpart[j],Point(Line[2],Line[3]),(float)mean_ht);
	  int typei = LineType(Line);
	  if(!((P1R1 == 1 && P2R1 == 1) || (P1R1 == 0 && P2R1 == 0)))
	  {
	    for(int m=0;m<contours_otherpart.size();m++)
	    {
	      if(ValidOtherParts[m] && m!=j)
	      {
		int P1R2 = NewPointRectangleTest(boundRect_otherpart[m],Point(Line[0],Line[1]),(float)mean_ht);
		int P2R2 = NewPointRectangleTest(boundRect_otherpart[m],Point(Line[2],Line[3]),(float)mean_ht);
		{
		  if(!((P1R2 == 1 && P2R2 == 1) || (P1R2 == 0 && P2R2 == 0)))
		  {		    
		    if(P1R1==1 && P2R2==1)// P1 in R1 and P2 in R2
		    {
		      Vec4i line1,line2;
		      if(typei == 1) // horizontal
		      {
			line1[1] = Line[1]; line1[3] = Line[1];
			line2[1] = Line[3]; line2[3] = Line[3];
			if(Line[2]>Line[0])
			{
			  line1[0] = Line[0];
			  line2[2] = Line[2];
			  line1[2] = boundRect_otherpart[j].x + boundRect_otherpart[j].width;
			  line2[0] = boundRect_otherpart[m].x;
			}
			else
			{
			  line1[0] = Line[0];
			  line2[2] = boundRect_otherpart[m].x + boundRect_otherpart[m].width;
			  line1[2] = boundRect_otherpart[j].x;
			  line2[0] = Line[2];
			}
		      }
		      if(typei == 2) // vertical
		      {
			line1[0] = Line[0]; line1[2] = Line[0];
			line2[0] = Line[2]; line2[2] = Line[2];
			if(Line[3]>Line[1])
			{
			  line1[1] = Line[1];
			  line2[3] = Line[3];
			  line1[3] = boundRect_otherpart[j].y + boundRect_otherpart[j].height;
			  line2[1] = boundRect_otherpart[m].y;
			}
			else
			{
			  line1[1] = Line[1];
			  line2[3] = boundRect_otherpart[m].y + boundRect_otherpart[m].height;
			  line1[3] = boundRect_otherpart[j].y;
			  line2[1] = Line[3];
			}
		      }
		      NewJoinedLine[k].push_back(line1);
		      NewJoinedLine[k].push_back(line2);
		      splitl = true;
		    }
		    if(P1R2==1 && P2R1==1)// P1 in R2 and P2 in R1
		    {
		      Vec4i line1,line2;
		      if(typei == 1) // horizontal
		      {
			line1[1] = Line[1]; line1[3] = Line[1];
			line2[1] = Line[3]; line2[3] = Line[3];
			if(Line[2]>Line[0])
			{
			  line1[0] = Line[0];
			  line2[2] = Line[2];
			  line1[2] = boundRect_otherpart[m].x + boundRect_otherpart[m].width;
			  line2[0] = boundRect_otherpart[j].x;
			}
			else
			{
			  line1[0] = Line[0];
			  line2[2] = boundRect_otherpart[j].x + boundRect_otherpart[j].width;
			  line1[2] = boundRect_otherpart[m].x;
			  line2[0] = Line[2];
			}
		      }
		      if(typei == 2) // vertical
		      {
			line1[0] = Line[0]; line1[2] = Line[0];
			line2[0] = Line[2]; line2[2] = Line[2];
			if(Line[3]>Line[1])
			{
			  line1[1] = Line[1];
			  line2[3] = Line[3];
			  line1[3] = boundRect_otherpart[m].y + boundRect_otherpart[m].height;
			  line2[1] = boundRect_otherpart[j].y;
			}
			else
			{
			  line1[1] = Line[1];
			  line2[3] = boundRect_otherpart[j].y + boundRect_otherpart[j].height;
			  line1[3] = boundRect_otherpart[m].y;
			  line2[1] = Line[3];
			}
		      }
		      if(typei == 1)
		      {
			if((max(line1[0],line1[2]) - min(line1[0],line1[2]))>2*mean_ht)
			  NewJoinedLine[k].push_back(line1);
			if((max(line2[0],line2[2]) - min(line2[0],line2[2]))>2*mean_ht)
			  NewJoinedLine[k].push_back(line2);
		      }
		      if(typei == 2)
		      {
			if((max(line1[1],line1[3]) - min(line1[1],line1[3]))>2*mean_ht)
			  NewJoinedLine[k].push_back(line1);
			if((max(line2[1],line2[3]) - min(line2[1],line2[3]))>2*mean_ht)
			  NewJoinedLine[k].push_back(line2);
		      }
		      splitl = true;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      if(!splitl)
	NewJoinedLine[k].push_back(JoinedLine[k][i]);
    }
  }
  
  
  for(size_t k = 0; k < NewJoinedLine.size(); k++)
  {
    Scalar color = Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) );
    for(int i=0;i<NewJoinedLine[k].size();i++)
    {
     
      line( LineImage, Point(NewJoinedLine[k][i][0], NewJoinedLine[k][i][1]),
	     Point(NewJoinedLine[k][i][2], NewJoinedLine[k][i][3]), color, 3, 8 );
    }
  }
  
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"NewJoinedLine.png");
  imwrite(name, LineImage);
  LineImage.release();
  
  JoinedLine.clear();
  JoinedLine = NewJoinedLine;
  NewJoinedLine.clear();
  
 /*********************************Working For Table Detection ************************************/ 
  
  vector<Rect> TempTableArea;
  vector<Vec4i> TableLines;
  SmoothedImg.copyTo(LineImage);
  cvtColor(LineImage,LineImage,CV_GRAY2BGR);
  for(size_t k = 0; k < JoinedLine.size(); k++)
  {
    for(int i=0;i<JoinedLine[k].size();i++)
    {
      int typei = LineType(JoinedLine[k][i]);
      if(typei == 1)
      {
	//cout << "ForLine " << k << "_" << i << "Working Lowerpart " << endl;
	vector<int> BP = TableHProjProfileLow(SmoothedImg,JoinedLine[k][i][1],min(JoinedLine[k][i][0],JoinedLine[k][i][2]),max(JoinedLine[k][i][0],JoinedLine[k][i][2]),mean_ht);
	//cout << "number of rows obtained is " << BP.size() << endl;
	if(BP.size()>0)
	{
	  TableLines.push_back(JoinedLine[k][i]);
	  line( LineImage, Point(JoinedLine[k][i][0], JoinedLine[k][i][1]),
	     Point(JoinedLine[k][i][2], JoinedLine[k][i][3]), Scalar(0,255,0), 3, 8 );
	//  cout << "ForLine " << k << "_" << i << "Working upperPartpart" << endl;
	  vector<int> TP = TableHProjProfileTop(SmoothedImg,JoinedLine[k][i][1],min(JoinedLine[k][i][0],JoinedLine[k][i][2]),max(JoinedLine[k][i][0],JoinedLine[k][i][2]),mean_ht);
	  int y_start = min(JoinedLine[k][i][1],JoinedLine[k][i][3])-TP.size();
	  int y_end = min(JoinedLine[k][i][1],JoinedLine[k][i][3])+BP.size();
	//  cout << "BP size=" << BP.size() << " TH size=" << TP.size() << endl;
	  rectangle(LineImage,Point(min(JoinedLine[k][i][0],JoinedLine[k][i][2]),y_start),Point(max(JoinedLine[k][i][0],JoinedLine[k][i][2]),y_end),Scalar(120,120,0));
	  Rect R(Point(min(JoinedLine[k][i][0],JoinedLine[k][i][2]),y_start),Point(max(JoinedLine[k][i][0],JoinedLine[k][i][2]),y_end));
	  TempTableArea.push_back(R);
	  //exit(0);
	//  vector<int> VP = TableVerticalProjection(SmoothedImg,y_start,y_end,min(JoinedLine[k][i][0],JoinedLine[k][i][2]),max(JoinedLine[k][i][0],JoinedLine[k][i][2]),mean_ht);
	}
      }
    }
  }
  vector<bool> TableFlag(TempTableArea.size(),true);
  
  /*src.copyTo(LineImage);
  for(int i=0;i<TempTableArea.size();i++)
  {
    if(TableFlag[i])
    {
      rectangle(LineImage,TempTableArea[i],Scalar(120,12,79));
    }
  }
  */
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"ProbableTableRegions_Lines.png");
  imwrite(name, LineImage);
  LineImage.release();
  
  
  
   // Combining two very near-by same size table (0.85)
  
  for(int i=0;i<TempTableArea.size();i++)
  {
    if(TableFlag[i])
    {
      for(int j=0;j<TempTableArea.size();j++)
      {
	if(TableFlag[j] && i!=j)
	{
	  int lenC = min(TempTableArea[i].x+TempTableArea[i].width,TempTableArea[j].x+TempTableArea[j].width) - max(TempTableArea[i].x,TempTableArea[j].x);
	  int lenU = max(TempTableArea[i].x+TempTableArea[i].width,TempTableArea[j].x+TempTableArea[j].width) - min(TempTableArea[i].x,TempTableArea[j].x);
	  float part_line = (lenC*1.0)/lenU;
	  float dist = GetHDistBetween2Rect(TempTableArea[i],TempTableArea[j]) * 1.0;
	  if(dist <= mean_ht && part_line>0.85)
	  {
	    Rect TmpR = GetUnionOf2Rect(TempTableArea[i],TempTableArea[j]);
	    TempTableArea[i] = TmpR;
	    TableFlag[j] = false;
	   // cout << "hello Rect " << j << " is removed and combined with " << i << endl;
	  }
	}
      }
    }
  }
  
  
  
  src.copyTo(LineImage);
  for(int i=0;i<TempTableArea.size();i++)
  {
    if(TableFlag[i])
    {
      rectangle(LineImage,TempTableArea[i],Scalar(120,12,79));
      
    }
  }
  
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"ProbableTableRegionsNewA.png");
  imwrite(name, LineImage);
  LineImage.release();
  
  
  
  
  for(int i=0;i<TempTableArea.size();i++)
  {
    if(TableFlag[i])
    {
      float area_sum = 0.0;
      float Common_pix_count = 0.0;
      for(int j=0;j<TempTableArea.size();j++)
      {
	if(i!=j && TableFlag[j])
	{
	  int val = FindOverlappingRectangles(TempTableArea[i],TempTableArea[j]);
	  float pixel_cnt_i = PixelsInaRect(SmoothedImg,TempTableArea[i])*1.0;
	  float pixel_cnt_j = PixelsInaRect(SmoothedImg,TempTableArea[j])*1.0;
	  if(val!=0)
	  {
	    if(val==1 || val ==2)
	    {
	      if(val==1)
	      {
		if((pixel_cnt_i/pixel_cnt_j)>0.75)
		{
		  if((TempTableArea[i].area()*1.0)/(TempTableArea[j].area()*1.0) > 0.65)
		    TableFlag[i] = false;
		  else
		    TableFlag[j] = false;
		}
		else
		  TableFlag[i] = false;
	      }
	      else
	      {
		if((pixel_cnt_j/pixel_cnt_i)>0.75)
		{
		  if((TempTableArea[j].area()*1.0)/(TempTableArea[i].area()*1.0) > 0.65)
		    TableFlag[j] = false;
		  else
		    TableFlag[i] = false;
		}
		else
		  TableFlag[j] = false;
	      }
	    }
	    else
	    {
	      Rect UR,IR;
	      UR = GetUnionOf2Rect(TempTableArea[i],TempTableArea[j]);
	      IR = GetInterSectionOf2Rect(TempTableArea[i],TempTableArea[j]);
	      area_sum = area_sum + (IR.area()*1.0);
	      Common_pix_count = Common_pix_count + (PixelsInaRect(SmoothedImg,IR)*1.0);
	      if(HLineRectangleTest(TempTableArea[i],TableLines[j],mean_ht)==1)
	      {
		  float h_per,v_per;
		  h_per = (min((TempTableArea[i].x+TempTableArea[i].width),max(TableLines[j][0],TableLines[j][2])) - max(TempTableArea[i].x,min(TableLines[j][0],TableLines[j][2])))*1.0;
		  h_per = h_per/(max(TempTableArea[i].width,(max(TableLines[j][0],TableLines[j][2])-min(TableLines[j][0],TableLines[j][2]))));
		  v_per = (min((TempTableArea[i].y+TempTableArea[i].height),(TempTableArea[j].y+TempTableArea[j].height)) - max(TempTableArea[i].y,TempTableArea[j].y))*1.0;
		  v_per =  v_per/max(TempTableArea[i].height,TempTableArea[j].height);
		  if(v_per>=0.8)
		  {
		    if(pixel_cnt_j/pixel_cnt_i>=0.85)
		    {
		      TableFlag[i] = false;
		    }
		  }
	      }
	      else
	      {
		  float h_per,v_per;
		  h_per = (min((TempTableArea[i].x+TempTableArea[i].width),max(TableLines[j][0],TableLines[j][2])) - max(TempTableArea[i].x,min(TableLines[j][0],TableLines[j][2])))*1.0;
		  h_per = h_per/(max(TempTableArea[i].width,(max(TableLines[j][0],TableLines[j][2])-min(TableLines[j][0],TableLines[j][2]))));
		  v_per = (min((TempTableArea[i].y+TempTableArea[i].height),(TempTableArea[j].y+TempTableArea[j].height)) - max(TempTableArea[i].y,TempTableArea[j].y))*1.0;
		  v_per =  v_per/max(TempTableArea[i].height,TempTableArea[j].height);
		  float hdist_L = (max(TempTableArea[i].x,TempTableArea[j].x) - min(TempTableArea[i].x,TempTableArea[j].x))*1.0;
		  if(hdist_L<3*mean_ht && v_per >=0.8)
		  {
		    if(pixel_cnt_j/pixel_cnt_i>=0.85)
		    {
		      TableFlag[i] = false;
		    }
		  }
	      }
	    }
	  }
	  
	  /*
	  if(val == 1 || val ==2)
	  {
	    if(val == 1)
	      TableFlag[i] = false;
	    else
	      TableFlag[j] = false;
	  }
	  else if(val == 3)
	  {
	      
	      Point2i TL1 = TempTableArea[i].tl(); Point2i TL2 = TempTableArea[j].tl();
	      Point2i BR1 = TempTableArea[i].br(); Point2i BR2 = TempTableArea[j].br();
	      Point2i TL,BR;
	      TL.x = max(TL1.x,TL2.x); TL.y = max(TL1.y,TL2.y);
	      BR.x = min(BR1.x,BR2.x); BR.y = min(BR1.y,BR2.y);
	      Rect IR(TL,BR);
	      TL.x = min(TL1.x,TL2.x); TL.y = min(TL1.y,TL2.y);
	      BR.x = max(BR1.x,BR2.x); BR.y = max(BR1.y,BR2.y);
	      Rect UR(TL,BR);
	      float ratio = (IR.area()*1.0)/min(TempTableArea[i].area(),TempTableArea[j].area());
	      if(ratio > 0.8)
	      {
		  TempTableArea[i] = UR;
		  TableFlag[j] = false;
	      }
	      else
	      {
		//Need to update	
		if(HLineRectangleTest(TempTableArea[i],TableLines[j],mean_ht)==1)
		{
		  float h_per,v_per;
		  h_per = (min((TempTableArea[i].x+TempTableArea[i].width),max(TableLines[j][0],TableLines[j][2])) - max(TempTableArea[i].x,min(TableLines[j][0],TableLines[j][2])))*1.0;
		  h_per = h_per/(max(TempTableArea[i].width,(max(TableLines[j][0],TableLines[j][2])-min(TableLines[j][0],TableLines[j][2]))));
		  v_per = (min((TempTableArea[i].y+TempTableArea[i].height),(TempTableArea[j].y+TempTableArea[j].height)) - max(TempTableArea[i].y,TempTableArea[j].y))*1.0;
		  v_per =  v_per/max(TempTableArea[i].height,TempTableArea[j].height);
		  if(v_per>=0.8)
		  {
		    TempTableArea[i] = UR;
		    TableFlag[j] = false;
		  }
		}
		else
		{
		  float h_per,v_per;
		  h_per = (min((TempTableArea[i].x+TempTableArea[i].width),max(TableLines[j][0],TableLines[j][2])) - max(TempTableArea[i].x,min(TableLines[j][0],TableLines[j][2])))*1.0;
		  h_per = h_per/(max(TempTableArea[i].width,(max(TableLines[j][0],TableLines[j][2])-min(TableLines[j][0],TableLines[j][2]))));
		  v_per = (min((TempTableArea[i].y+TempTableArea[i].height),(TempTableArea[j].y+TempTableArea[j].height)) - max(TempTableArea[i].y,TempTableArea[j].y))*1.0;
		  v_per =  v_per/max(TempTableArea[i].height,TempTableArea[j].height);
		  float hdist_L = (max(TempTableArea[i].x,TempTableArea[j].x) - min(TempTableArea[i].x,TempTableArea[j].x))*1.0;
		  if(hdist_L<3*mean_ht && v_per >=0.8)
		  {
		    TempTableArea[i] = UR;
		    TableFlag[j] = false;
		  }
		}
		
	     }	      
	  }
	  */
	  
	} // End of If in J
      } // End of j
      float TotalPixelCnt = PixelsInaRect(SmoothedImg,TempTableArea[i])*1.0;
      float part_own_area = (TempTableArea[i].area() - area_sum)/TempTableArea[i].area();
      float per_own_pixcnt = TotalPixelCnt - Common_pix_count;
      per_own_pixcnt = per_own_pixcnt/TotalPixelCnt;
      //cout << "For Rect " << i << " area overlapped is " << part_own_area << endl;
      if(part_own_area<0.5 && per_own_pixcnt < 0.15)
	TableFlag[i] = false;
      //if(no_overlap>0)
	//TableFlag[i] = false;
    }
  }
  
  /*
  src.copyTo(LineImage);
  for(int i=0;i<TempTableArea.size();i++)
  {
    if(TableFlag[i])
    {
      rectangle(LineImage,TempTableArea[i],Scalar(120,12,79));
    }
  }
  
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"ProbableTableRegionsNewB.png");
  imwrite(name, LineImage);
  LineImage.release();
  
  
  
  for(int i=0;i<TempTableArea.size();i++)
  {
    if(TableFlag[i])
    {
      float area_sum = 0.0;
      float Common_pix_count = 0.0;
      for(int j=0;j<TempTableArea.size();j++)
      {
	if(i!=j && TableFlag[j])
	{
	  int val = FindOverlappingRectangles(TempTableArea[i],TempTableArea[j]);
	  if(val == 3)
	  {
	    Point2i TL1 = TempTableArea[i].tl(); Point2i TL2 = TempTableArea[j].tl();
	    Point2i BR1 = TempTableArea[i].br(); Point2i BR2 = TempTableArea[j].br();
	    Point2i TL,BR;
	    TL.x = max(TL1.x,TL2.x); TL.y = max(TL1.y,TL2.y);
	    BR.x = min(BR1.x,BR2.x); BR.y = min(BR1.y,BR2.y);
	    Rect IR(TL,BR);
	    area_sum = area_sum + (IR.area()*1.0);
	    Common_pix_count = Common_pix_count + (PixelsInaRect(SmoothedImg,IR)*1.0);
	  }
	}
      }
      float TotalPixelCnt = PixelsInaRect(SmoothedImg,TempTableArea[i])*1.0;
      float part_own_area = (TempTableArea[i].area() - area_sum)/TempTableArea[i].area();
      float per_own_pixcnt = TotalPixelCnt - Common_pix_count;
      per_own_pixcnt = per_own_pixcnt/TotalPixelCnt;
      //cout << "For Rect " << i << " area overlapped is " << part_own_area << endl;
      if(part_own_area<0.5 && per_own_pixcnt < 0.15)
	TableFlag[i] = false;
    }
  }
  */
  src.copyTo(LineImage);
  for(int i=0;i<TempTableArea.size();i++)
  {
    if(TableFlag[i])
    {
      rectangle(LineImage,TempTableArea[i],Scalar(120,12,79));
    }
  }
  
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"ProbableTableRegions_BasedonLines.png");
  imwrite(name, LineImage);
  LineImage.release();
  
  
  
  
  
  //vector<bool> ValidOtherParts(contours_otherpart.size(),true);
  vector<Rect> OtherPart_Rect;
  vector<bool> worked(contours_otherpart.size(),true);
  for( int j = 0; j < contours_otherpart.size(); j++ )
  {
    if(ValidOtherParts[j] && worked[j])
    {
      for( int i = 0; i < contours_otherpart.size(); i++ )
      {
	if(!ValidOtherParts[i])
	  worked[i] = false;
	if(ValidOtherParts[i] && i!=j && worked[i])
	{	  
	  int val = FindOverlappingRectangles(boundRect_otherpart[j],boundRect_otherpart[i]);
	  if(val==3)
	  {
	    Point2i TL1 =  boundRect_otherpart[j].tl(); Point2i TL2 = boundRect_otherpart[i].tl();
	    Point2i BR1 =  boundRect_otherpart[j].br(); Point2i BR2 = boundRect_otherpart[i].br();
	    Point2i TL,BR;
	    TL.x = min(TL1.x,TL2.x); TL.y = min(TL1.y,TL2.y);
	    BR.x = max(BR1.x,BR2.x); BR.y = max(BR1.y,BR2.y);
	    Rect UR(TL,BR);
	    OtherPart_Rect.push_back(UR);
	    worked[j] = false;
	    worked[i] = false;
	  }
	}
      }
    }
  }
  for( int j = 0; j < contours_otherpart.size(); j++ )
  {
    if(worked[j])
    {
      OtherPart_Rect.push_back(boundRect_otherpart[j]);
      worked[j] = false;
    }
  }
 /* for(int i=0;i<TempTableArea.size();i++)
  {
    if(TableFlag[i])
    {
      for( int j = 0; j < contours_otherpart.size(); j++ )
      {	
	if(hierarchy_otherpart[j][3] == -1 && ValidOtherParts[j])
	{
	  int val1 = FindOverlappingRectangles(TempTableArea[i],boundRect_otherpart[j]);
	  if(val1!=0)
	  {
	    for(int k = 0;k<contours_otherpart.size();k++)
	    {
	      if(j!=k)
	      {
		if(hierarchy_otherpart[k][3] == -1 && ValidOtherParts[k])
		{
		  int val2 = FindOverlappingRectangles(TempTableArea[i],boundRect_otherpart[k]);
		  if(val2!=0)
		  {
		    //Combine Rect j and k
		    Point2i TL1 =  boundRect_otherpart[j].tl(); Point2i TL2 = boundRect_otherpart[k].tl();
		    Point2i BR1 =  boundRect_otherpart[j].br(); Point2i BR2 = boundRect_otherpart[k].br();
		    Point2i TL,BR;
		    TL.x = min(TL1.x,TL2.x); TL.y = min(TL1.y,TL2.y);
		    BR.x = max(BR1.x,BR2.x); BR.y = max(BR1.y,BR2.y);
		    Rect UR(TL,BR);
		    OtherPart_Rect.push_back(UR);
		    ValidOtherParts[j] = false;
		    ValidOtherParts[k] = false;
		  }
		}
		else
		  ValidOtherParts[k] = false;
	      }
	    }
	  }
	}
	else
	  ValidOtherParts[j] = false;
      }
    }
  }
  
  for( int j = 0; j < contours_otherpart.size(); j++ )
  {
    if(ValidOtherParts[j])
      OtherPart_Rect.push_back(boundRect_otherpart[j]);
  }*/
  
  src.copyTo(LineImage);
  for(int i=0;i<OtherPart_Rect.size();i++)
  {
    rectangle(LineImage,OtherPart_Rect[i],Scalar(120,12,79));
  }
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"bounded_Otherpart_new.png");
  imwrite(name, LineImage);
  LineImage.release();
  
  for(int i=0;i<TempTableArea.size();i++)
  {
    if(TableFlag[i])
    {
      cout << "Working for Rect " << i << endl;
      cout << TempTableArea[i].y << "," << TempTableArea[i].x << "::" << TempTableArea[i].height << "," << TempTableArea[i].width << endl;
      for( int j = 0; j < OtherPart_Rect.size(); j++ )
      {
	cout << "Working for LC " << j << endl;
	cout << OtherPart_Rect[j].y  << "," << OtherPart_Rect[j].x  << "::" << OtherPart_Rect[j].height  << "," << OtherPart_Rect[j].width << endl;
	int val = FindOverlappingRectangles(TempTableArea[i],OtherPart_Rect[j]);
	cout << "RectangleTest Value= " << val << endl;
	if(val!=0)
	  cout << "RectangleTest Value= " << val << endl;
	float ratio;
	if(val == 1 || val ==2)
	{	    
	  if(val == 1)
	  {
	    ratio = (TempTableArea[i].area()*1.0)/OtherPart_Rect[j].area();
	    cout << "raio= " << ratio << endl;
	    if(ratio >= 0.6)
	    {
	      TempTableArea[i] = OtherPart_Rect[j];
	    }
	    else
	    {
	      TableFlag[i] = false;
	      cout << "Due to ratio condition rejecting the rect " << i  << endl; 
	    }
	  }
	  else
	  {
	    ratio = (OtherPart_Rect[j].area()*1.0)/TempTableArea[i].area();
	    cout << "raio= " << ratio << endl;
	  }
	}
	else if(val == 3)
	{
	  Rect IR = GetInterSectionOf2Rect(TempTableArea[i],OtherPart_Rect[j]);
	  Rect UR = GetUnionOf2Rect(TempTableArea[i],OtherPart_Rect[j]);	  
	  ratio = (IR.area()*1.0)/UR.area();
	  cout << "raio= " << ratio << endl;
	  if(OtherPart_Rect[j].height>=3*mean_ht)
	  {
	    cout << "I am in with height " << OtherPart_Rect[j].height << endl;
	    if(ratio>=0.55)
	    {
	      TempTableArea[i] = UR;
	    }
	    else
	    {
	      float new_ratio = (IR.area()*1.0)/OtherPart_Rect[j].area();
	      float h_per,v_per;
	      h_per = min(TempTableArea[i].x+TempTableArea[i].width,OtherPart_Rect[j].x+OtherPart_Rect[j].width) - max(TempTableArea[i].x,OtherPart_Rect[j].x); 
	      h_per = h_per/OtherPart_Rect[j].width;
	      v_per = min(TempTableArea[i].y+TempTableArea[i].height,OtherPart_Rect[j].y+OtherPart_Rect[j].height) - max(TempTableArea[i].y,OtherPart_Rect[j].y);
	      v_per = v_per/OtherPart_Rect[j].height;
	      if(h_per>0.7 && v_per>0.7)
	      {
		TempTableArea.push_back(OtherPart_Rect[j]);
		TableFlag.push_back(true);
	      }
	      else
		cout << "In overlapping Condition rejecting the rect " << i  << endl;
	      TableFlag[i] = false;
	    }
	    
	  }
	}
      }
    }
  }
  
  src.copyTo(LineImage);
  for(int i=0;i<TempTableArea.size();i++)
  {
    if(TableFlag[i])
    {
      rectangle(LineImage,TempTableArea[i],Scalar(120,12,79));
    }
  }
  
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"ProbableTableRegions_B.png");
  imwrite(name, LineImage);
  LineImage.release();
  
  for(int i=0;i<TempTableArea.size();i++)
  {
    if(TableFlag[i])
    {
      for( int j = 0; j < OtherPart_Rect.size(); j++ )
      {
	int val = FindOverlappingRectangles(TempTableArea[i],OtherPart_Rect[j]);
	if(val == 0)
	{
	  int dist = GetHDistBetween2Rect(TempTableArea[i],OtherPart_Rect[j]);
	  if(dist<=1.5*mean_ht && OtherPart_Rect[j].height < 1.5*mean_ht)
	  {
	    Rect UR = GetUnionOf2Rect(TempTableArea[i],OtherPart_Rect[j]);
	    TempTableArea[i] = UR;
	  }
	}
      }
    }
  }
  
  src.copyTo(LineImage);
  for(int i=0;i<TempTableArea.size();i++)
  {
    if(TableFlag[i])
    {
      rectangle(LineImage,TempTableArea[i],Scalar(120,12,79));
    }
  }
  
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"ProbableTableRegions_FineTuned.png");
  imwrite(name, LineImage);
  LineImage.release();
  
  
  hgap = (int) ((mean_ht*0.75)+0.5);
  SmoothedImg = horizontal_gapfilling(NewBinaryDst,hgap);
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"Word_inBinary.png");
  imwrite(name, SmoothedImg);
  
  
  vector<Rect> NewBR = FindBoundRects(SmoothedImg,Rect(0,0,SmoothedImg.cols,SmoothedImg.rows));
  Mat NewSImage = Mat(SmoothedImg.rows,SmoothedImg.cols,SmoothedImg.type(),Scalar(255));
  
  for(int i=0;i<NewBR.size();i++)
  {
    Rect rr = NewBR[i];
    if(rr.height<3*mean_ht)
    {
      for(int r=rr.y;r<rr.y+rr.height;r++)
      {
	for(int c=rr.x;c<rr.x+rr.width;c++)
	  NewSImage.at<uchar>(r,c) = SmoothedImg.at<uchar>(r,c);
      }
    }
  }
  
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"Word_inBinary_Filtered.png");
  imwrite(name, NewSImage);
  
  NewSImage.copyTo(SmoothedImg);
  
   // Combining two very near-by same size table (0.85)
  
  for(int i=0;i<TempTableArea.size();i++)
  {
    if(TableFlag[i])
    {
      for(int j=0;j<TempTableArea.size();j++)
      {
	if(TableFlag[j] && i!=j)
	{
	  int lenC = min(TempTableArea[i].x+TempTableArea[i].width,TempTableArea[j].x+TempTableArea[j].width) - max(TempTableArea[i].x,TempTableArea[j].x);
	  int lenU = max(TempTableArea[i].x+TempTableArea[i].width,TempTableArea[j].x+TempTableArea[j].width) - min(TempTableArea[i].x,TempTableArea[j].x);
	  float part_line = (lenC*1.0)/lenU;
	  float dist = GetHDistBetween2Rect(TempTableArea[i],TempTableArea[j]) * 1.0;
	  if(dist <= mean_ht && part_line>0.85)
	  {
	    Rect TmpR = GetUnionOf2Rect(TempTableArea[i],TempTableArea[j]);
	    TempTableArea[i] = TmpR;
	    TableFlag[j] = false;
	    cout << "hello Rect " << j << " is removed and combined with " << i << endl;
	  }
	}
      }
    }
  }
  
  
  
  src.copyTo(LineImage);
  for(int i=0;i<TempTableArea.size();i++)
  {
    if(TableFlag[i])
    {
      rectangle(LineImage,TempTableArea[i],Scalar(120,12,79));
      
    }
  }
  
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"ProbableTableRegionsNew.png");
  imwrite(name, LineImage);
  LineImage.release();
  
  
  
  for(int i=0;i<TempTableArea.size();i++)
  {
    if(TableFlag[i])
    {
      Rect R = TempTableArea[i];
    //  cout << "Working for PT " << i << endl;
      vector<int> VP = TableVerticalProjection(SmoothedImg,R.y,R.y+R.height,R.x,R.x+R.width,mean_ht);
      int nzcnt = 0;
      int zcnt = 0;
      for(int j=0;j<VP.size();j++)
      {
	if(VP[j]<0)
	  nzcnt++;
	else
	  zcnt++;
      }
    //  cout << "There are " << nzcnt << " non-zero elemens" << endl;
      if(nzcnt<2)
	TableFlag[i] = false;
    }
  }
  
  
   src.copyTo(LineImage);
  for(int i=0;i<TempTableArea.size();i++)
  {
    if(TableFlag[i])
    {
      cout << "print Table " << i << " with x,y=" << TempTableArea[i].x << "," << TempTableArea[i].y << endl;
      rectangle(LineImage,TempTableArea[i],Scalar(120,12,79));
      
    }
  }
  
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"ProbableTableRegions.png");
  imwrite(name, LineImage);
  LineImage.release();
  
  
  
  
  
 
  
  FILE *fp;
  fp = fopen("submission.xml","a+");
  fprintf(fp,"<document filename=\"%s\">\n",argv[1]);
  src.copyTo(LineImage);
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder("LabeledImage",argv[1]);
  src.copyTo(LineImage);
  Mat TableImage = imread(name,1);
  for(int i=0;i<TempTableArea.size();i++)
  {
    if(TableFlag[i])
    {
      cout << "For Table " << i << " with x,y=" << TempTableArea[i].x << "," << TempTableArea[i].y << endl;
      if(GetStructureofaTable(SmoothedImg,TempTableArea[i],(int) mean_ht+0.5))
      {
	//FineTuneTheTableRegion(NewBinaryDst,TempTableArea[i],(float) mean_ht);
	if(FineTuneTheTableRegion(SmoothedImg, TempTableArea[i], OtherPart_Rect, JoinedLine, (int) (mean_ht+0.5)))
	{
	  rectangle(LineImage,TempTableArea[i],Scalar(120,12,79),3);
	  rectangle(TableImage,TempTableArea[i],Scalar(100,170,80),3);
	  Rect R = TempTableArea[i];
	  fprintf(fp,"\t<tableRegion prob=\"1.0\">\n");
	  fprintf(fp,"\t\t<Coords points=\"%d,%d %d,%d %d,%d %d,%d\"/>\n",R.x,R.y,R.x+R.width,R.y,R.x,R.y+R.height,R.x+R.width,R.y+R.height);
	  fprintf(fp,"\t</tableRegion>\n");
	}
	else
	  TableFlag[i] = false;
      }
      else
	TableFlag[i] = false;
    }
  }
  
  
  /*
  hgap = (int) ((mean_ht*1.5)+0.5);
  SmoothedImg = horizontal_gapfilling(NewBinaryDst,hgap);
  
  vector<Rect> NewTables = FindTablesWithNoHLines(SmoothedImg,TempTableArea,(int) (mean_ht+0.5));
  
  for(int i=0;i<NewTables.size();i++)
  {
	rectangle(LineImage,NewTables[i],Scalar(120,12,79),3);
	rectangle(TableImage,NewTables[i],Scalar(100,170,80),3);
	Rect R = NewTables[i];
	fprintf(fp,"\t<tableRegion prob=\"1.0\">\n");
	fprintf(fp,"\t\t<Coords points=\"%d,%d %d,%d %d,%d %d,%d\"/>\n",R.x,R.y,R.x+R.width,R.y,R.x,R.y+R.height,R.x+R.width,R.y+R.height);
	fprintf(fp,"\t</tableRegion>\n");
  }
  
  */
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"TableRegions.png");
  imwrite(name, LineImage);
  LineImage.release();
  
  makedir("TableImages");
  string timgname(substring);
  timgname = timgname + ".png";
  tempname = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  strcpy(tempname,timgname.c_str());
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder("TableImages",tempname);
  imwrite(name, TableImage);
  TableImage.release();
  
  
  
  
  /************** Checking for Proper Multicolumn **************/
  
  cout << "Not Working on Figure and Formula" << endl;
  /*
  vector<bool> valid_OtherPart_Rect(OtherPart_Rect.size(),true);
  
  for(int j = 0; j < OtherPart_Rect.size(); j++)
  {
    if(valid_OtherPart_Rect[j])
    {
      for(int i = 0; i < OtherPart_Rect.size(); i++)
      {
	if(i!=j && valid_OtherPart_Rect[i])
	{
	  int val = FindOverlappingRectangles(OtherPart_Rect[j],OtherPart_Rect[i]);
	  if(val!=0)
	  {
	    if(val == 1)
	    {
	      valid_OtherPart_Rect[j] = false;
	    }
	    else if(val == 2)
	      valid_OtherPart_Rect[i] = false;
	    else
	    {
		Point2i TL1 = OtherPart_Rect[i].tl(); Point2i TL2 = OtherPart_Rect[j].tl();
		Point2i BR1 = OtherPart_Rect[i].br(); Point2i BR2 = OtherPart_Rect[j].br();
		Point2i TL,BR;
		TL.x = min(TL1.x,TL2.x); TL.y = min(TL1.y,TL2.y);
		BR.x = max(BR1.x,BR2.x); BR.y = max(BR1.y,BR2.y);
		Rect UR(TL,BR);
		if(OtherPart_Rect[i].area()>OtherPart_Rect[j].area())
		{
		  OtherPart_Rect[i] = UR;
		  valid_OtherPart_Rect[j] = false;
		}
		else
		{
		  OtherPart_Rect[j] = UR;
		  valid_OtherPart_Rect[i] = false;
		}
	    }
	  } 
	}
      }
    }
  }
  
  src.copyTo(LineImage);
  
  
  for(int j = 0; j < OtherPart_Rect.size(); j++)
  {
    if(valid_OtherPart_Rect[j])
    {
      bool flag = true;
      for(int i=0;i<TempTableArea.size();i++)
      {
	if(TableFlag[i])
	{
	  int val = FindOverlappingRectangles(TempTableArea[i],OtherPart_Rect[j]);
	  if(val!=0)
	    flag = false;
	}
      }
      if(flag && OtherPart_Rect[j].height>=1.5*mean_ht)
      {
	rectangle(LineImage,OtherPart_Rect[j],Scalar(120,12,79));
	Rect R = OtherPart_Rect[j];
	fprintf(fp,"\t<figureRegion prob=\"1.0\">\n");
	fprintf(fp,"\t\t<Coords points=\"%d,%d %d,%d %d,%d %d,%d\"/>\n",R.x,R.y,R.x+R.width,R.y,R.x,R.y+R.height,R.x+R.width,R.y+R.height);
	fprintf(fp,"\t</figureRegion>\n");
      }
    }
  }
  
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"ProbableGraphicsRegions.png");
  imwrite(name, LineImage);
  LineImage.release();
  */
  
  fprintf(fp,"</document>\n");
  fclose(fp);
  /*
  for(int i=0;i<TempTableArea.size();i++)
  {
    if(TableFlag[i])
    {
      cout << "Working for Rect " << i << endl;
      cout << TempTableArea[i].y << "," << TempTableArea[i].x << "::" << TempTableArea[i].height << "," << TempTableArea[i].width << endl;
      for( int j = 0; j < contours_otherpart.size(); j++ )
      {
	if(hierarchy_otherpart[j][3] == -1)
	{
	  cout << "Working for LC " << j << endl;
	  cout << boundRect_otherpart[j].y  << "," << boundRect_otherpart[j].x  << "::" << boundRect_otherpart[j].height  << "," << boundRect_otherpart[j].width << endl;
	  int val = FindOverlappingRectangles(TempTableArea[i],boundRect_otherpart[j]);
	  cout << "RectangleTest Value= " << val << endl;
	  if(val!=0)
	    cout << "RectangleTest Value= " << val << endl;
	  float ratio;
	  if(val == 1 || val ==2)
	  {	    
	    if(val == 1)
	    {
	      ratio = (TempTableArea[i].area()*1.0)/boundRect_otherpart[j].area();
	      cout << "raio= " << ratio << endl;
	      if(ratio >= 0.75)
	      {
		TempTableArea[i] = boundRect_otherpart[j];
	      }
	      else
		TableFlag[i] = false;
	    }
	    else
	    {
	      ratio = (boundRect_otherpart[j].area()*1.0)/TempTableArea[i].area();
	      cout << "raio= " << ratio << endl;
	    }
	  }
	  else if(val == 3)
	  {
	    Point2i TL1 = TempTableArea[i].tl(); Point2i TL2 = boundRect_otherpart[j].tl();
	    Point2i BR1 = TempTableArea[i].br(); Point2i BR2 = boundRect_otherpart[j].br();
	    Point2i TL,BR;
	    TL.x = max(TL1.x,TL2.x); TL.y = max(TL1.y,TL2.y);
	    BR.x = min(BR1.x,BR2.x); BR.y = min(BR1.y,BR2.y);
	    Rect IR(TL,BR);
	    TL.x = min(TL1.x,TL2.x); TL.y = min(TL1.y,TL2.y);
	    BR.x = max(BR1.x,BR2.x); BR.y = max(BR1.y,BR2.y);
	    Rect UR(TL,BR);
	    ratio = (IR.area()*1.0)/UR.area();
	    cout << "raio= " << ratio << endl;
	    if(boundRect_otherpart[j].height>=1.5*mean_ht)
	    {
	      cout << "I am in with height " << boundRect_otherpart[j].height << endl;
	      if(ratio>=0.55)
	      {
		TempTableArea[i] = UR;
	      }
	      else
		TableFlag[i] = false;
	    }
	  }
	}
      }
    }
  }
  
  src.copyTo(LineImage);
  for(int i=0;i<TempTableArea.size();i++)
  {
    if(TableFlag[i])
    {
      rectangle(LineImage,TempTableArea[i],Scalar(120,12,79));
    }
  }
  
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"ProbableTableRegions_FineTuned.png");
  imwrite(name, LineImage);
  LineImage.release();
  */
  
  
  
  
  
 /* 
  src.copyTo(LineImage);
  
  vector<vector<Vec4i> > ProbableParas;
  
  for(size_t k = 0; k < JoinedLine.size(); k++)
  {  
    vector<bool> LineFlag(JoinedLine[k].size(),true);
    
    for(int i=0;i<JoinedLine[k].size();i++)
    {
      int typei = LineType(JoinedLine[k][i]);
      
      if(LineFlag[i] && typei == 1)
      {
	int max_y,min_y;
	Vec4i Line_1,Line_2;
	Line_1 = JoinedLine[k][i];
	Line_2 = JoinedLine[k][i];
	cout << "y_i="  << min(JoinedLine[k][i][1],JoinedLine[k][i][3]) << endl;
	max_y = min(JoinedLine[k][i][1],JoinedLine[k][i][3]);
	min_y = max_y;
	int cnt_sameline = 0;
	LineFlag[i] = false;
	Scalar color = Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) );
      // line( LineImage, Point(JoinedLine[i][0], JoinedLine[i][1]),
	  //   Point(JoinedLine[i][2], JoinedLine[i][3]), color, 3, 8 );
	for(int j=0;j<JoinedLine[k].size();j++)
	{
	  int typej = LineType(JoinedLine[k][j]);
	  if(typej == typei && LineFlag[j])
	  {
	    double frac = partCoveredByTwoParallelLine(JoinedLine[k][i],JoinedLine[k][j]);
	    printf("Frac between Line %d and %d is %lf\n",i,j,frac);
	    if(frac > 0.9)
	    {
	      if(max_y < min(JoinedLine[k][j][1],JoinedLine[k][j][3]))
	      {
		max_y = min(JoinedLine[k][j][1],JoinedLine[k][j][3]);
		Line_2 = JoinedLine[k][j];
		//cout << "in Max y are from orig " << JoinedLine[k][j][1] << JoinedLine[k][j][3] << endl;
		//cout << "in Max y are from copied " << Line_2[1] << Line_2[3] << endl;
	      }
	      if(min_y > min(JoinedLine[k][j][1],JoinedLine[k][j][3]))
	      {
		min_y = min(JoinedLine[k][j][1],JoinedLine[k][j][3]);
		Line_1 = JoinedLine[k][j];
		//cout << "in Min y are from orig " << JoinedLine[k][j][1] << JoinedLine[k][j][3] << endl;
		//cout << "in Min y are from copied " << Line_1[1] << Line_1[3] << endl;
	      }
	     // cout << "y_j="  << min(JoinedLine[k][j][1],JoinedLine[k][j][3]) << endl;
	     // cout << "max=" << max_y << "min=" << min_y << endl;
	     // cout << "Assigning " <<"ymax=" << Line_2[1] << Line_2[3] << "ymin=" << Line_1[1] << Line_1[3] << endl;
	      cnt_sameline++;
	      LineFlag[j] = false;
	      line( LineImage, Point(JoinedLine[k][i][0], JoinedLine[k][i][1]),
	     Point(JoinedLine[k][i][2], JoinedLine[k][i][3]), color, 3, 8 );
	     line( LineImage, Point(JoinedLine[k][j][0], JoinedLine[k][j][1]),
	     Point(JoinedLine[k][j][2], JoinedLine[k][j][3]), color, 3, 8 );
	     //cout << "Number of overlapping Lines are " << cnt_sameline << endl;
	    }
	  }
	}
	vector<Vec4i> reg;
	if(cnt_sameline==0)
	{
	  line( LineImage, Point(Line1[0], Line1[1]),
	     Point(Line1[2], Line1[3]), color, 3, 8 );
	  reg.push_back(Line_1);
	 // cout << "max_new=" << min(Line_1[1],Line_1[3]) << "min_new=" << min(Line_1[1],Line_1[3]) << endl;
	}
	else
	{
	  reg.push_back(Line_1);
	  reg.push_back(Line_2);
	  line( LineImage, Point(Line1[0], Line1[1]),
	     Point(Line1[2], Line1[3]), color, 3, 8 );
	  line( LineImage, Point(Line2[0], Line2[1]),
	     Point(Line2[2], Line2[3]), color, 3, 8 );
	// cout << "ymax=" << Line_2[1] << Line_2[3] << "ymin=" << Line_1[1] << Line_1[3] << endl;
	// cout << "max_new=" << min(Line_2[1],Line_2[3]) << "min_new=" << min(Line_1[1],Line_1[3]) << endl;
	}
	ProbableParas.push_back(reg);  
      }
      
    }
    LineFlag.clear();
  }
  
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"LineCheck.png");
  imwrite(name, LineImage);
  LineImage.release();
  
  src.copyTo(LineImage);
  
  for(int i=0;i<ProbableParas.size();i++)
  {
    Scalar color = Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) );
   // cout << "Number of Lines are " << ProbableParas[i].size() << endl;
    for(int j=0;j<ProbableParas[i].size();j++)
    {
      Vec4i Line = ProbableParas[i][j];
      line( LineImage, Point(Line[0], Line[1]),
	     Point(Line[2], Line[3]), color, 3, 8 );
    }
  }
  
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"ProbableParaLines.png");
  imwrite(name, LineImage);
  LineImage.release();
  
  
  Mat HGImage = horizontal_gapfilling(NewBinaryDst,8);
  Mat VGImage = vertical_gapfilling(HGImage,5);
  
  name = (char *) malloc ( 2001 * sizeof(char));
  if(name == NULL)
  {
    printf("Memory can not be allocated\n");
    exit(0);
  }
  name = CreateNameIntoFolder(substring,"GapfilledImage.png");
  imwrite(name, VGImage);
  
  contours.clear();
  contours_poly.clear();
  boundRect.clear();
  boundRotatedRect.clear();
  hierarchy.clear();
  
  if(contours.empty())
    printf("Contor is empty\n");
  
  VGImage.copyTo(TempImage);
  VGImage.release();
  
  TempImage = FindImageInverse(TempImage);
  
  /// Find contours
  findContours( TempImage, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, Point(0, 0) );	
  
  contours_poly.resize(contours.size());
  boundRect.resize(contours.size());
  boundRotatedRect.resize(contours.size());
  
  for( int j = 0; j < contours.size(); j++ )
  { 
    approxPolyDP( Mat(contours[j]), contours_poly[j], 3, true );
    boundRect[j] = boundingRect( Mat(contours_poly[j]) );
    boundRotatedRect[j] = minAreaRect( Mat(contours_poly[j]) );
  }
  
  */
  
  
  return 0;
}
