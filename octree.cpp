//#include "stdafx.h"
#include "StdAfx.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "octree.h"
using namespace std;
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

Vec3 makeVec3( double x, double y, double z)
{
  //Vec3 v3 = (Vec3) malloc(3 * sizeof(double));//改成NEW就没有bug了，不知道为什么
  Vec3 v3=new double;
  v3[0] = x; v3[1] = y; v3[2] = z;
  return v3;
}

Vec3 copyVec3( Vec3 src )
{
 //Vec3 v3 = (Vec3) malloc(3 * sizeof(double));
  Vec3 v3=new double;
  v3[0] = src[0]; v3[1] = src[1]; v3[2] = src[2];
  return v3;
}

/* ------------------------------------------------------------------------ */

Octree* Octree::make_octree( Vec3 min, Vec3 max )
{
  //Octree* o = (Octree*) malloc(sizeof(Octree));
   Octree* o = new Octree;
  o->min = copyVec3(min);
  o->max = copyVec3(max);
  o->children = 0;
  o->at_max_depth = 0;
  /*
  printf("creating octree %.3lf,%.3lf,%.3lf ... %.3lf,%.3lf,%.3lf\n", 
	 o->min[2], o->min[1], o->min[0],
	 o->max[2], o->max[1], o->max[0] );
  */
  return o;
}
void Octree::subpoint( Octree* o,int oc,Vec3 min, Vec3 max)
{
	pvex_nor *m_p1,*m_p2;
	POSITION pos1,pos2;
	for(pos1=o->vex.GetHeadPosition(),pos2=o->normal.GetHeadPosition();pos1!=NULL;)
	{
		//pos1=o->vex.FindIndex(i);//pos2=o->normal.FindIndex(i);
		m_p1=(pvex_nor*)o->vex.GetNext(pos1);m_p2=(pvex_nor*)o->normal.GetNext(pos2);
		if((m_p1->x>min[0]&&m_p1->x<max[0])&&(m_p1->y>min[1]&&m_p1->y<max[1])
			&&(m_p1->z>min[2]&&m_p1->z<max[2]))
		{
			o->children[oc]->vex.AddHead(new pvex_nor(m_p1->x,m_p1->y,m_p1->z,m_p1->ID));
			o->children[oc]->normal.AddHead(new pvex_nor(m_p2->x,m_p2->y,m_p2->z,m_p2->ID));
		}
	}

}
void Octree::split_octree( Octree* o )
{

  double oc_min[3];
  double oc_max[3];
  
  Vec3 mid = makeVec3( (o->min[0] + o->max[0]) * 0.5, 
		       (o->min[1] + o->max[1]) * 0.5, 
		       (o->min[2] + o->max[2]) * 0.5 );
  int xp, yp, zp;
  int oc = 0;

  //o->children = (Octree**) malloc( 8 * sizeof(Octree*));
  o->children = new Octree*;
  for(zp=0; zp < 2; zp++)
  {
    if(zp == 0)
    { 
      oc_min[2] = o->min[2];
      oc_max[2] = mid[2];
    }
    else
    {
      oc_min[2] = mid[2];
      oc_max[2] = o->max[2];
    }
    
    for(yp=0; yp < 2; yp++)
    {
      if(yp == 0)
      { 
	oc_min[1] = o->min[1];
	oc_max[1] = mid[1];
      }
      else
      {
	oc_min[1] = mid[1];
	oc_max[1] = o->max[1];
      }
      
      for(xp=0; xp < 2; xp++)
      {
	if(xp == 0)
	{ 
	  oc_min[0] = o->min[0];
	  oc_max[0] = mid[0];
	}
	else
	{
	  oc_min[0] = mid[0];
	  oc_max[0] = o->max[0];
	}
	
	o->children[ (zp*4) + (yp*2) + xp ] = make_octree( oc_min, oc_max );
    subpoint( o,(zp*4) + (yp*2) + xp,oc_min,oc_max);
      }
    }
  }

}

/* ------------------------------------------------------------------------ */


int Octree::recursively_evaluate_octree( int min_depth, int max_depth, int current_depth, Octree* o )
{
  int deepest_child = current_depth;
  int oc=0;
  int cd=0;
  o->density=evaluate1_point(o);
  current_depth++;
  o->not_fully_divided = (char) octree_needs_to_be_split( o );
  if( current_depth <= max_depth )
  {
    if(( current_depth <= min_depth) || ( o->not_fully_divided ))
    {
	      split_octree( o );
      for(oc = 0; oc < 8; oc++)
      {
	     /*Vex.RemoveAll();
         for( pos = vex[oc].GetHeadPosition(); pos != NULL; )
		 {
                m_p=(pvex_nor*)vex[oc].GetNext( pos );
	        	Vex.AddHead(new pvex_nor(m_p->x,m_p->y,m_p->z));
		 }*/
         //if(deepest_child==current_depth||deepest_child==0)
		 //{
	     cd = recursively_evaluate_octree( min_depth, max_depth, current_depth, o->children[ oc ] );
		 //}
    
     	if(cd > deepest_child)
	    deepest_child = cd;
      }
    }
  }
  else
  {
    o->at_max_depth = 1;
  }
  
  return deepest_child;
}

/* ------------------------------------------------------------------------ */

int Octree::subdivide_octree( int min_depth, int max_depth, Octree* o )
{
  return recursively_evaluate_octree(min_depth, max_depth, 0, o );
}

double demo1( Vec3 pos )
{
  /* demo 1: the surface is a sphere of radius 0.75 centered at 0,0,0 
     
     function returns 1.0 if point inside sphere, 0.0 otherwise 
  */
  
  double dist_sq = (pos[0] * pos[0]) + (pos[1] * pos[1]) + (pos[2] * pos[2]);
  
  return ( dist_sq < 0.5625 ) ? 1.0 : 0.0;
}

double demo2( Vec3 pos )
{
  /* demo 2: the surface is two spheres, 
           A: radius 0.25 centered at -.25,-.5,0 
       and B: radius 0.5  centered at -0.5,0,0 

     function returns 1.0 if point inside sphere A, 2.0 for sphere B, 0.0 for neither
  */

  double dist_sq_a = ((pos[0]+.25) * (pos[0]+.25)) + ((pos[1]+.5) * (pos[1]+.5)) + (pos[2] * pos[2]);
  double dist_sq_b = ((pos[0]+.8) * (pos[0]+.8)) + (pos[1] * pos[1]) + (pos[2] * pos[2]);

  if( dist_sq_a <= .0625 )
    return 1.0;


  if( dist_sq_b <= .25 )
    return 2.0;

  return 0.0;
}

double demo3( Vec3 pos )
{
  /* demo 3: the surface is tiny sphere, radius 0.1 centered at -.5,.5,0 

     function returns 1.0 if point inside sphere A, 0.0 otherwise
  */


  double dist_sq = ((pos[0]+.5) * (pos[0]+.5)) + ((pos[1]-.5) * (pos[1]-.5)) + (pos[2] * pos[2]);

  return ( dist_sq < 0.01 ) ? 1.0 : 0.0;

}

double demo4( Vec3 pos )
{
  /* demo 4: wavey surface
     
     function returns 1.0 if point 'near' surface , 0.0 otherwise
     
  */
  double surface_height = sin( (pos[0] * 3.0) ) * cos ( (pos[1] * 3.0) );
  
  double distance_sq = (pos[2] - surface_height) * (pos[2] - surface_height);

  return ( distance_sq < 0.01 ) ? 1.0 : 0.0;
}

double demo5( Vec3 pos )
{
  /* demo 5: hemisphere, center 0,0,0 radius 0.5, cut by plane at z=0

   */

  double abs_dist_sq = ((pos[0]) * (pos[0])) + ((pos[1]) * (pos[1])) + (pos[2] * pos[2]);
  
  double surf_dist_sq = abs_dist_sq - 0.5625;
  if(surf_dist_sq  < 0)
    surf_dist_sq = -surf_dist_sq;

  if( (pos[2] > 0) && (surf_dist_sq < 0.1 ))
    return 1.0;
  else
    return .0;
}

double demo6( Vec3 pos )
{
  /* demo 6: another wavey surface
     
     function returns 1.0 if point 'near' surface , 0.0 otherwise
     
  */
  double surface_height = sin( (pos[0] * 2.0) ) + sin ( (pos[1] * 2.0) );
  
  double distance_sq = (pos[2] - surface_height) * (pos[2] - surface_height);

  return ( distance_sq < 0.01 ) ? 1.0 : 0.0;
}


double demo7( Vec3 pos )
{
  /* demo 7: a cylinder
     
     function returns 1.0 if point 'near' surface , 0.0 otherwise
     
  */

  double disc_dist_sq = ((pos[1]) * (pos[1])) + ((pos[2] * pos[2]));

  double surf_dist_sq =  disc_dist_sq - 0.5625;
   if(surf_dist_sq  < 0)
    surf_dist_sq = -surf_dist_sq;
  
  if( ( pos[0] > -0.75 ) && ( pos[0] < 0.75 ) )
  {
    if(  surf_dist_sq < 0.1 )
      return 1.0;
  }

  return 0.0;
}


double demo8( Vec3 pos )
{
  /* demo8 : mandlebrot */
  int max_iters = 500;

  int it_count = 0;

  double cr, ci, zr, zi, new_zr, new_zi;
  int inside = 1;

  /* only exists near the plane z=0 */
  
  if((pos[2] > 0.001) || (pos[2] < -0.001))
    return .0;

  zr = cr = (pos[0] * 2.0);
  zi = ci = (pos[1] * 2.0);

  while((it_count < max_iters) && (inside))
  {
    /* z = z^2 + c */

    it_count++;

    new_zr = (zr*zr) - (zi*zi) + cr;
    new_zi = (2*zr*zi) + ci;

    zr = new_zr;
    zi = new_zi;

    inside = (zr * zr) + (zi * zi) < (2.0 * 20);
  }

  return (it_count == max_iters) ? 1.0 : 0.0;
  

}
double demo9( Vec3 pos, Octree *o)
{
	pvex_nor *m_p1,*m_p2;
	POSITION po;
	double dis;
	double zx,zy,zz,tempx,tempy,tempz,temp;
   
    for(int i=0;i<o->vex.GetCount();i++)
	{   
		po=o->vex.FindIndex(i);
		m_p1=(pvex_nor *)o->vex.GetAt(po);
		m_p2=(pvex_nor *)o->normal.GetAt(po);
		tempx=pos[0]-m_p1->x;tempy=pos[1]-m_p1->y;tempz=pos[2]-m_p1->z;
		temp=tempx*m_p2->x+tempy*m_p2->y+tempz*m_p2->z;
		zx=m_p1->x-temp*m_p2->x;zy=m_p1->y-temp*m_p2->y;zz=m_p1->z-temp*m_p2->z;
		dis=sqrt((pos[0]-zx)*(pos[0]-zx)+(pos[1]-zy)*(pos[1]-zy)+(pos[2]-zz)*(pos[2]-zz));
		//dis=sqrt((pos[0]-m_p1->x)*(pos[0]-m_p1->x)+(pos[1]-m_p1->y)*(pos[1]-m_p1->y)
			//+(pos[2]-m_p1->y)*(pos[2]-m_p1->y));
		if(dis<0.1) goto loop;
	}
loop:
	return (dis<0.1) ? 1.0 : 0.0;
}

/* ------------------------------------------------------------------------ */

double Octree::evaluate_point( Vec3 pos,Octree *o )
{
  return demo9( pos,o );
}
int Octree::evaluate1_point( Octree *o )
{
	int tCnt=o->vex.GetCount()+1;
	return(tCnt);

}

/* ------------------------------------------------------------------------ */

/* returns 1 if the octree should be split, 0 otherwise */

/*
 this implementation checks whether all 8 corner values are the same

 (if any corner 1..7 is different to corner 0 then the function returns 1)
*/

int Octree::octree_needs_to_be_split( Octree* o )
{
  /*int i;
  double v = o->value[0];

  for(i=1; i < 8; i++)
    if( o->value[i] != v)
	     return 1;

  /* if we got here, then all corners have the same value */

  //return 0;
  if(o->density>20) return 1;
  else return 0;
}
void Octree::isoface(Octree* o)
{   
 //	char bf[2];
	////CPtrList vexlist,norlist;
	//float v1,v2,v3;
	//ifstream file("ww9.obj");
	//if(!o->vex.IsEmpty())
	//	o->vex.RemoveAll();
	//if(!o->normal.IsEmpty())
	//	o->normal.RemoveAll(); 
 //  	if(!file.fail())
	//{
	//	while(!file.eof())
	//	{
	//		file>>bf>>v1>>v2>>v3;
	//		if(bf[0]=='v'&&bf[1]==0)
	//			o->vex.AddTail(new pvex_nor(v1,v2,v3));
	//		if(bf[0]=='v'&&bf[1]=='n')
	//			o->normal.AddTail(new pvex_nor(v1,v2,v3));
	//	}
	//	//o->vex=vexlist;
	//}

 //   file.close();
}

/*
   Linearly interpolate the position where an isosurface cuts
   an edge between two vertices, each with their own scalar value
*/
Vec3 VertexInterp(double isolevel,Vec3 p1,Vec3 p2,double valp1,double valp2)
{
   double mu;
   Vec3 p= makeVec3(0,0,0);

   if (fabs(isolevel-valp1) < 0.00001)
      return(p1);
   if (fabs(isolevel-valp2) < 0.00001)
      return(p2);
   if (fabs(valp1-valp2) < 0.00001)
      return(p1);
   mu = (isolevel - valp1) / (valp2 - valp1);
   p[0] = p1[0] + mu * (p2[0] - p1[0]);
   p[1] = p1[1] + mu * (p2[1] - p1[1]);
   p[2] = p1[2] + mu * (p2[2] - p1[2]);

   return(p);
}
/*
   Given a grid cell and an isolevel, calculate the triangular
   facets required to represent the isosurface through the cell.
   Return the number of triangular facets, the array "triangles"
   will be loaded up with the vertices at most 5 triangular facets.
	0 will be returned if the grid cell is either totally above
   of totally below the isolevel.
*/
int Polygonise(GRIDCELL grid,double isolevel,TRIANGLE *triangles)
{
   int i,ntriang;
   int cubeindex;
   Vec3 vertlist[12];

   int edgeTable[256]={
   0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
   0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
   0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
   0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
   0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
   0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
   0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
   0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
   0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
   0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
   0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
   0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
   0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
   0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
   0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
   0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
   0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
   0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
   0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
   0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
   0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
   0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
   0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
   0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
   0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
   0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
   0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
   0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
   0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
   0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
   0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
   0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };
int triTable[256][16] =
{{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
{3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
{3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
{3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
{9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
{9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
{2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
{8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
{9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
{4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
{3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
{1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
{4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
{4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
{5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
{2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
{9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
{0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
{2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
{10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
{5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
{5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
{9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
{0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
{1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
{10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
{8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
{2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
{7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
{2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
{11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
{5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
{11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
{11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
{1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
{9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
{5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
{2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
{5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
{6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
{3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
{6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
{5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
{1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
{10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
{6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
{8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
{7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
{3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
{5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
{0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
{9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
{8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
{5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
{0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
{6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
{10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
{10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
{8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
{1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
{0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
{10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
{3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
{6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
{9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
{8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
{3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
{6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
{0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
{10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
{10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
{2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
{7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
{7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
{2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
{1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
{11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
{8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
{0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
{7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
{10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
{2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
{6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
{7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
{2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
{1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
{10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
{10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
{0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
{7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
{6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
{8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
{9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
{6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
{4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
{10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
{8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
{0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
{1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
{8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
{10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
{4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
{10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
{5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
{11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
{9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
{6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
{7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
{3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
{7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
{3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
{6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
{9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
{1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
{4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
{7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
{6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
{3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
{0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
{6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
{0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
{11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
{6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
{5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
{9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
{1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
{1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
{10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
{0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
{5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
{10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
{11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
{9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
{7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
{2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
{8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
{9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
{9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
{1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
{9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
{9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
{5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
{0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
{10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
{2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
{0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
{0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
{9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
{5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
{3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
{5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
{8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
{0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
{9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
{1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
{3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
{4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
{9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
{11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
{11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
{2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
{9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
{3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
{1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
{4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
{3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
{0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
{9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
{1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

   /*
      Determine the index into the edge table which
      tells us which vertices are inside of the surface
   */
   cubeindex = 0;
   if (grid.val[0] < isolevel) cubeindex |= 1;
   if (grid.val[1] < isolevel) cubeindex |= 2;
   if (grid.val[2] < isolevel) cubeindex |= 4;
   if (grid.val[3] < isolevel) cubeindex |= 8;
   if (grid.val[4] < isolevel) cubeindex |= 16;
   if (grid.val[5] < isolevel) cubeindex |= 32;
   if (grid.val[6] < isolevel) cubeindex |= 64;
   if (grid.val[7] < isolevel) cubeindex |= 128;

   /* Cube is entirely in/out of the surface */
   if (edgeTable[cubeindex] == 0)
      return(0);

   /* Find the vertices where the surface intersects the cube */
   if (edgeTable[cubeindex] & 1)
      vertlist[0] =
         VertexInterp(isolevel,grid.p[0],grid.p[1],grid.val[0],grid.val[1]);
   if (edgeTable[cubeindex] & 2)
      vertlist[1] =
         VertexInterp(isolevel,grid.p[1],grid.p[2],grid.val[1],grid.val[2]);
   if (edgeTable[cubeindex] & 4)
      vertlist[2] =
         VertexInterp(isolevel,grid.p[2],grid.p[3],grid.val[2],grid.val[3]);
   if (edgeTable[cubeindex] & 8)
      vertlist[3] =
         VertexInterp(isolevel,grid.p[3],grid.p[0],grid.val[3],grid.val[0]);
   if (edgeTable[cubeindex] & 16)
      vertlist[4] =
         VertexInterp(isolevel,grid.p[4],grid.p[5],grid.val[4],grid.val[5]);
   if (edgeTable[cubeindex] & 32)
      vertlist[5] =
         VertexInterp(isolevel,grid.p[5],grid.p[6],grid.val[5],grid.val[6]);
   if (edgeTable[cubeindex] & 64)
      vertlist[6] =
         VertexInterp(isolevel,grid.p[6],grid.p[7],grid.val[6],grid.val[7]);
   if (edgeTable[cubeindex] & 128)
      vertlist[7] =
         VertexInterp(isolevel,grid.p[7],grid.p[4],grid.val[7],grid.val[4]);
   if (edgeTable[cubeindex] & 256)
      vertlist[8] =
         VertexInterp(isolevel,grid.p[0],grid.p[4],grid.val[0],grid.val[4]);
   if (edgeTable[cubeindex] & 512)
      vertlist[9] =
         VertexInterp(isolevel,grid.p[1],grid.p[5],grid.val[1],grid.val[5]);
   if (edgeTable[cubeindex] & 1024)
      vertlist[10] =
         VertexInterp(isolevel,grid.p[2],grid.p[6],grid.val[2],grid.val[6]);
   if (edgeTable[cubeindex] & 2048)
      vertlist[11] =
         VertexInterp(isolevel,grid.p[3],grid.p[7],grid.val[3],grid.val[7]);

   /* Create the triangle */
   ntriang = 0;
   for (i=0;triTable[cubeindex][i]!=-1;i+=3) {
      triangles[ntriang].p[0] = vertlist[triTable[cubeindex][i  ]];
      triangles[ntriang].p[1] = vertlist[triTable[cubeindex][i+1]];
      triangles[ntriang].p[2] = vertlist[triTable[cubeindex][i+2]];
      ntriang++;
   }

   return(ntriang);
}
/*marchingcube*/
double signal_distance(Vec3 pos, Octree *o)
{
	pvex_nor *m_p1,*m_p2;
	POSITION position1,position2;
	double dis,dismin=99999.0;
	double gv=0;
	double zx,zy,zz,tempx,tempy,tempz,temp;
   
    for(position1=o->vex.GetHeadPosition(),position2=o->normal.GetHeadPosition();position1!=NULL;)
	{   
		m_p1=(pvex_nor *)o->vex.GetNext(position1);
		m_p2=(pvex_nor *)o->normal.GetNext(position2);
		tempx=pos[0]-m_p1->x;tempy=pos[1]-m_p1->y;tempz=pos[2]-m_p1->z;
		temp=tempx*m_p2->x+tempy*m_p2->y+tempz*m_p2->z;
		zx=m_p1->x-temp*m_p2->x;zy=m_p1->y-temp*m_p2->y;zz=m_p1->z-temp*m_p2->z;
		dis=sqrt((pos[0]-zx)*(pos[0]-zx)+(pos[1]-zy)*(pos[1]-zy)+(pos[2]-zz)*(pos[2]-zz));
		if(dismin>dis) {dismin=dis;gv=temp;}
	}
	return (gv);

}
void marchingcube(int depth,Octree *o)
{   
	Vec3 point = makeVec3(0,0,0);
	int xp, yp, zp,c, l = 0;
	int nCount;
	double isolevel;
	GRIDCELL grid;
	TRIANGLE *triangles= new TRIANGLE;
	//TRIANGLE *triangles = (TRIANGLE*) malloc(sizeof(TRIANGLE));
    isolevel=0;
	if(o->children==0&&o->vex.GetCount()>5)
	{
	    for(zp=0; zp < 2; zp++)
		{
             point[2] = (zp == 0) ? o->min[2] : o->max[2];

            for(yp=0; yp < 2; yp++)
			{
				point[1] = (yp == 0) ? o->min[1] : o->max[1];

                for(xp=0; xp < 2; xp++)
				{ 
	                point[0] = (xp == 0) ? o->min[0] : o->max[0];
	                grid.p[l] = point;
		            grid.val[l]=signal_distance(point,o);
					l++;
				}
			}
		}	
        nCount=Polygonise(grid,isolevel,triangles);
	}
	if(o->children!=0)
	{
		for(c=0; c < 8; c++)
		    marchingcube( depth + 1, o->children[ c ] );
	  
    }
  
}

///////////////////////////////八叉树求交点
///////////////////////////////1.建立八叉树
void Octree::creat_octree(Octree*& pOctree,vec_PFACETTRI &vec_face,vector<POINT3D>&  vec_box)
{
	double p1[3],p2[3],dp;
	p1[0]=0;p1[1]=0;p1[2]=0;
	p2[0]=0;p2[1]=0;p2[2]=0;
	dp=0;
	dp=vec_box[1].x-vec_box[0].x;
	if ((vec_box[1].y-vec_box[0].y)>dp)
	{
		dp=vec_box[1].y-vec_box[0].y;
	}
	if ((vec_box[1].z-vec_box[0].z)>dp)
	{
		dp=vec_box[1].z-vec_box[0].z;
	}
	p1[0]=vec_box[0].x;p1[1]=vec_box[0].y;p1[2]=vec_box[0].z;
	p2[0]=vec_box[0].x+dp;p2[1]=vec_box[0].y+dp;p2[2]=vec_box[0].z+dp;
	pOctree=pOctree->make_octree(p1,p2);
	for (int i=0;i<vec_face.size();i++)
	{
		POINT3D temV=Ccalcubase::Center(vec_face[i]->pHEdge);
		pOctree->vex.AddTail(new pvex_nor(temV.x,temV.y,temV.z,i));
		pOctree->normal.AddTail(new pvex_nor(vec_face[i]->m_PFacetNorm->dx,vec_face[i]->m_PFacetNorm->dy,vec_face[i]->m_PFacetNorm->dz,i));
	}
	subdivide_octree(0,8,pOctree);
}
void Octree::delet_octree_1(Octree*& pOctree)
{
	queue<Octree*> Q;
	Q.push(pOctree);
	while(1)
	{
		int num=0;
		for (int i=0;i<8;i++)
		{
			if (Q.front()->children[i]->not_fully_divided==1)
			{
				num++;
				Q.push(Q.front()->children[i]);
				Q.pop();
				break;
			}
		}
		if (num==0)
		{
			for (int i=0;i<8;i++)
			{
				Octree* ptemp=Q.front()->children[i];
				if (ptemp!=NULL)
				{
					delete ptemp;
				}
				Q.front()->children[i]=NULL;
				ptemp=NULL;
			}
			Octree* ptemp=Q.front();
			delete ptemp;
			ptemp=NULL;
			break;
		}	
	}
}
void Octree::delet_octree(Octree*& pOctree)
{
	while(1)
	{
		if (pOctree->not_fully_divided==1)
		{
			delet_octree_1(pOctree);
		}
		if (pOctree->not_fully_divided==0)
		{
			delete pOctree;
			pOctree=NULL;
			break;
		}
	}
}
BOOL Octree::get_face(Octree* pOctree,VECTOR3D LineNor,POINT3D orig,PFACETTRI faceA,vector<int>& ID_crossface)
{
	BOOL has=0;
	queue<Octree*> Q_o;
	vector<int> vec_crossface;
	POINT3D retu;
	Q_o.push(pOctree);
	while (!Q_o.empty())
	{
		for (int i=0;i<8;i++)
		{
			if (Q_o.front()->not_fully_divided==1)
				if (cros_BigCube(Q_o.front()->children[i],LineNor,orig)==1)
				{
					Q_o.push(Q_o.front()->children[i]);
				}
		}
		if (Q_o.front()->not_fully_divided==0)
		{
			has=TRUE;
			POSITION po;
			pvex_nor *m_p;
			for (int i=0;i<Q_o.front()->vex.GetCount();i++)
			{
				po=Q_o.front()->vex.FindIndex(i);
				m_p=(pvex_nor *)Q_o.front()->vex.GetAt(po);
				ID_crossface.push_back(m_p->ID);
			}
		}
		Q_o.pop();
	}
	return has;
}
PFACETTRI Octree::get_face_new(Octree* &pOctree,VECTOR3D LineNor,POINT3D orig,PFACETTRI faceA,vec_PFACETTRI& m_vecPFacetTri)//回溯法遍历
{
	queue<Octree*> Q_o;
	Q_o.push(pOctree);
	while (!Q_o.empty())
	{
		for (int i=0;i<8;i++)
		{
			if (Q_o.front()->not_fully_divided==1)
				if (cros_BigCube(Q_o.front()->children[i],LineNor,orig)==1)
				{
					Q_o.push(Q_o.front()->children[i]);
				}
		}
		if (Q_o.front()->not_fully_divided==0)
		{
			POSITION po;
			pvex_nor *m_p;
			for (int i=0;i<Q_o.front()->vex.GetCount();i++)
			{
				po=Q_o.front()->vex.FindIndex(i);
				m_p=(pvex_nor *)Q_o.front()->vex.GetAt(po);
				PFACETTRI theface=m_vecPFacetTri[m_p->ID];
				POINT3D V[3];
				V[0]=(POINT3D)(*theface->m_PVerts[0]);
				V[1]=(POINT3D)(*theface->m_PVerts[1]);
				V[2]=(POINT3D)(*theface->m_PVerts[2]);
				if (Ccalcubase::IntersectTriangle(orig,LineNor,V[0],V[1],V[2])&&theface!=faceA)
				{
					return theface;
				}
			}
		}
		Q_o.pop();
	}
	return NULL;
}
PFACETTRI Octree::get_face_new1(Octree* &pOctree,VECTOR3D LineNor,POINT3D orig,PFACETTRI faceA,vec_PFACETTRI& m_vecPFacetTri)//回溯法遍历
{
	queue<Octree*> Q_o;
	Q_o.push(pOctree);
	while (!Q_o.empty())
	{
		for (int i=0;i<8;i++)
		{
			if (Q_o.front()->not_fully_divided==1)
				if (cros_BigCube(Q_o.front()->children[i],LineNor,orig)==1)
				{
					Q_o.push(Q_o.front()->children[i]);
				}
		}
		if (Q_o.front()->not_fully_divided==0)
		{
			POSITION po;
			pvex_nor *m_p;
			for (int i=0;i<Q_o.front()->vex.GetCount();i++)
			{
				po=Q_o.front()->vex.FindIndex(i);
				m_p=(pvex_nor *)Q_o.front()->vex.GetAt(po);
				PFACETTRI theface=m_vecPFacetTri[m_p->ID];
				POINT3D V[3];
				V[0]=(POINT3D)(*theface->m_PVerts[0]);
				V[1]=(POINT3D)(*theface->m_PVerts[1]);
				V[2]=(POINT3D)(*theface->m_PVerts[2]);
				VECTOR3D v1;
				double v1n=0;
				v1=(VECTOR3D)(*theface->m_PFacetNorm);
				v1n=v1|LineNor;
				if (Ccalcubase::IntersectTriangle(orig,LineNor,V[0],V[1],V[2])&&theface->bStatus==0&&v1n>0.8)
				{
					return theface;
				}
			}
		}
		Q_o.pop();
	}
	return NULL;
}
BOOL Octree::cros_BigCube(Octree* pOctree,VECTOR3D LineNor,POINT3D orig)
{
	POINT3D v[8];
     /*
	 父节点包围盒的八个顶点顺序
	     6---7
        /|  /|
       5---4 |
       | 3-|-2
       |/  |/
       0---1
       */
	v[0].x=pOctree->min[0];v[0].y=pOctree->min[1];v[0].z=pOctree->min[2];
	v[1].x=pOctree->max[0];v[1].y=pOctree->min[1];v[1].z=pOctree->min[2];
	v[2].x=pOctree->max[0];v[2].y=pOctree->max[1];v[2].z=pOctree->min[2];
	v[3].x=pOctree->min[0];v[3].y=pOctree->max[1];v[3].z=pOctree->min[2];
	v[4].x=pOctree->max[0];v[4].y=pOctree->min[1];v[4].z=pOctree->max[2];
	v[5].x=pOctree->min[0];v[5].y=pOctree->min[1];v[5].z=pOctree->max[2];
	v[6].x=pOctree->min[0];v[6].y=pOctree->max[1];v[6].z=pOctree->max[2];
	v[7].x=pOctree->max[0];v[7].y=pOctree->max[1];v[7].z=pOctree->max[2];
	if (Ccalcubase::IntersectTriangle(orig,LineNor,v[0],v[1],v[2]))
		return 1;
	else if (Ccalcubase::IntersectTriangle(orig,LineNor,v[0],v[2],v[3]))
	    return 1;
	else if (Ccalcubase::IntersectTriangle(orig,LineNor,v[0],v[1],v[4]))
		return 1;
	else if (Ccalcubase::IntersectTriangle(orig,LineNor,v[0],v[4],v[5]))
	    return 1;
	else if (Ccalcubase::IntersectTriangle(orig,LineNor,v[1],v[2],v[4]))
	    return 1;
	else if (Ccalcubase::IntersectTriangle(orig,LineNor,v[2],v[4],v[7]))
	    return 1;
	else if (Ccalcubase::IntersectTriangle(orig,LineNor,v[0],v[3],v[5]))
	    return 1;
	else if (Ccalcubase::IntersectTriangle(orig,LineNor,v[3],v[5],v[6]))
	    return 1;
	else if (Ccalcubase::IntersectTriangle(orig,LineNor,v[2],v[3],v[7]))
	    return 1;
	else if (Ccalcubase::IntersectTriangle(orig,LineNor,v[3],v[6],v[7]))
	    return 1;
	else if (Ccalcubase::IntersectTriangle(orig,LineNor,v[4],v[5],v[6]))
		return 1;
	else if (Ccalcubase::IntersectTriangle(orig,LineNor,v[4],v[6],v[7]))
		return 1;
	else
		return FALSE;
}
