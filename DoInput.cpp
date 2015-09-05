// DoInput.cpp : 实现文件
//

#include "stdafx.h"
#include "WZHViewer.h"
#include "DoInput.h"
#include "afxdialogex.h"


// DoInput 对话框

IMPLEMENT_DYNAMIC(DoInput, CDialogEx)

DoInput::DoInput(CWnd* pParent /*=NULL*/)
	: CDialogEx(DoInput::IDD, pParent)
	, theLIMIT(2.6)
{

}

DoInput::~DoInput()
{
}

void DoInput::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDIT2, theLIMIT);
}


BEGIN_MESSAGE_MAP(DoInput, CDialogEx)
END_MESSAGE_MAP()


// DoInput 消息处理程序
