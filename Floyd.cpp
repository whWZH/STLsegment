#include "StdAfx.h"
#include "Floyd.h"
CFloyd::CFloyd(void)
{

}
CFloyd::~CFloyd(void)
{
}
 void CFloyd::init(int n)
 {
	 floyd_Bigmap.resize(n);
	 for (int i=0;i<n;i++)
	 {
		 floyd_Bigmap[i].resize(n);
	 }
	const double INF = 90000;
    for(int i=0;i<n;i++)
	{
        for(int j=0;j<n;j++)
		{
           if (i==j)
           {
			   floyd_Bigmap[i][j]=0;
           }
		   else
		   {
			   floyd_Bigmap[i][j]=INF;
		   }
		}
	}
 }
 void CFloyd::floyd(int n)
 {
	 floyd_dist.resize(n);
	 for (int i=0;i<n;i++)
	 {
		 floyd_dist[i].resize(n);
	 }
    for(int i=0;i<n;i++)
	{
	 for(int j=0;j<n;j++)
	 {
        floyd_dist[i][j]=floyd_Bigmap[i][j]/*,path[i][j]=0*/;
	 }
	}
      for(int k=0;k<n;k++)
	  {
        for(int i=0;i<n;i++)
		{
          for(int j=0;j<n;j++)
		  {
            if(floyd_dist[i][k]+floyd_dist[k][j]<floyd_dist[i][j])
			{
         	  floyd_dist[i][j]=floyd_dist[i][k]+floyd_dist[k][j]/*,path[i][j]=k*/;
			}
		  }
		}
	  }
 }
