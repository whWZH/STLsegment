// Bin2AscOfSTL.h : main header file for the BIN2ASCOFSTL application
//

#if !defined(AFX_BIN2ASCOFSTL_H__FD7CAE7C_E2D6_4EA7_A033_81E349930EB5__INCLUDED_)
#define AFX_BIN2ASCOFSTL_H__FD7CAE7C_E2D6_4EA7_A033_81E349930EB5__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"		// main symbols

/////////////////////////////////////////////////////////////////////////////
// CBin2AscOfSTLApp:
// See Bin2AscOfSTL.cpp for the implementation of this class
//

class CBin2AscOfSTLApp : public CWinApp
{
public:
	CBin2AscOfSTLApp();

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CBin2AscOfSTLApp)
	public:
	virtual BOOL InitInstance();
	//}}AFX_VIRTUAL

// Implementation

	//{{AFX_MSG(CBin2AscOfSTLApp)
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_BIN2ASCOFSTL_H__FD7CAE7C_E2D6_4EA7_A033_81E349930EB5__INCLUDED_)
