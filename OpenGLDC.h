#pragma once
#include "afx.h"

#include "gl/GL.h"
#include "gl/GLU.h"
#include "gl/glut.h"
#include "gl/glaux.h"

#include "GCamera.h"
#include "GeomBase.h"
#include "mesh.h"

#define BUFFER_LENGTH 64

class COpenGLDC :public CObject
{
public:
	COpenGLDC(HWND hWnd);
	virtual~COpenGLDC(void);

private:
	HWND   m_hWnd;       //关联窗口的句柄
	HGLRC  m_hRC;        //渲染场境句柄
	HDC    m_hDC;        //设备场境句柄

	COLORREF m_clrBk;       //窗口背景色
	COLORREF m_clr;         //非光照模式下的模型颜色
	COLORREF m_clrHighlight;//用于高亮显示时的模型颜色，如拾取一个物体时需要高亮度显示
	COLORREF m_clrMaterial; //材料的颜色

	BOOL     m_bShading;    //是否采用着色显示
	GLfloat  m_vecLight[3]; //光源方向

//    BOOL     m_bSelectionMode;            //当前是否是选择模式
//	  GLuint   m_selectBuff[BUFFER_LENGTH]; //选择缓存区

public:
	GCamera  m_Camera;    //照相机，用于取景操作

protected:
	void ClearBKground(); //清除背景颜色
	void OnShading();     //光照/非光照模式设置

public:
	BOOL InitDC();                    //初始化
	void GLResize(int cx,int cy);     //对应窗口尺寸变化
	void GLSetupRC();                 //设置渲染场境

	void Ready();                     //绘制前准备函数
	void Finish();                    //结束绘图函数

////////////////////对光照与颜色的操作函数//////////////////////

	void Shading(BOOL bShading);        //光照/非光照模式切换
	BOOL IsShading();                   //当前是否是着色模式

	void Lighting(BOOL bLighting);      //是否使用光源
	BOOL IsLighting();
	BOOL m_bSlelectionMode;
	GLuint m_selectBuff[BUFFER_LENGTH];

	void SetLightDirection(float dx,float dy,float dz);  //设置与获取光源方向
	void GetLightDirection(float& dx,float& dy,float& dz);

	void SetMaterialColor(COLORREF clr);            //设置与获取材料颜色
	void GetMaterialColor(COLORREF& clr);

	void SetBkColor(COLORREF rgb);                 //设置与获取背景颜色
	void GetBkColor(COLORREF& rgb);

	void SetColor(COLORREF rgb);                  //设置与获取非光照模式下的绘制颜色
	void GetColor(COLORREF& rgb);

	void SetHighlightColor(COLORREF clr);         //设置与获取高亮度显示的颜色
	void GetHighlightColor(COLORREF& clr);

	void Highlight(BOOL bLight = TRUE);           //高亮度/正常显示切换
	POINT3D BeginSelection(int xPos,int yPos);
	int EndSelection(UINT* items);
	BOOL isSelectionMode();
	/*void InitNames();
	void LoadName(UINT name);
	void PushName(UINT name);
	void PopName();
*/

////////////////////////基于OpenGl基本图元的高级绘图函数//////////////////////////////
	void DrawPoint(const POINT3D&);            //绘制一个空间点
	void DrawCoord();                          //绘制用户坐标系
	void DrawLine(const POINT3D& sp,const POINT3D& ep);//绘制一条直线
	void DrawLine2(const POINT3D& sp,const POINT3D& ep);
	void DrawTriChip(double n0,double n1,double n2,double v00,double v01,double v02,       //绘制一个三角面片
		             double v10,double v11,double v12,double v20,double v21,double v22);
	/////////////////////////////////
	/*static void RenderInformation(int txtA,int txtB);
	static  void DrawText(char* string, int flag=18);*/

};
