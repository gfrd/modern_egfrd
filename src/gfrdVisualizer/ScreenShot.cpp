#if defined(_MSC_VER)
#include <windows.h>

#pragma warning( push )
#pragma warning( disable : 4458)
#include <gdiplus.h>
#pragma warning( pop )

#pragma warning( push )
#pragma warning( disable : 4505)
#include <GL/freeglut.h>
#pragma warning( pop )

// --------------------------------------------------------------------------------------------------------------------------------

static ULONG_PTR gdiplusToken;

// --------------------------------------------------------------------------------------------------------------------------------

void screen_capture_startup()
{
   Gdiplus::GdiplusStartupInput gdiplusStartupInput;
   Gdiplus::GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, nullptr);
}

// --------------------------------------------------------------------------------------------------------------------------------

void screen_capture_shutdown()
{
   Gdiplus::GdiplusShutdown(gdiplusToken);
}

// --------------------------------------------------------------------------------------------------------------------------------

int GetEncoderClsid(const WCHAR* format, CLSID* pClsid)
{
   using namespace Gdiplus;
   UINT  num = 0;          // number of image encoders
   UINT  size = 0;         // size of the image encoder array in bytes

   GetImageEncodersSize(&num, &size);
   if (size == 0) return -1;

   ImageCodecInfo* pImageCodecInfo = static_cast<ImageCodecInfo*>(malloc(size));
   if (pImageCodecInfo == nullptr) return -1;

   GetImageEncoders(num, size, pImageCodecInfo);
   for (UINT j = 0; j < num; ++j)
   {
      if (wcscmp(pImageCodecInfo[j].MimeType, format) == 0)
      {
         *pClsid = pImageCodecInfo[j].Clsid;
         free(pImageCodecInfo);
         return j;
      }
   }

   free(pImageCodecInfo);
   return 0;
}

// --------------------------------------------------------------------------------------------------------------------------------

void screenshot(int number)
{
   int width = glutGet(GLUT_WINDOW_WIDTH);
   int height = glutGet(GLUT_WINDOW_HEIGHT);

   int stride = ((width * 24 + 31) / 32) * 4;
   int size = stride * height;
   unsigned char *data = new unsigned char[size];
   glReadPixels(0, 0, width, height, GL_BGR_EXT, GL_UNSIGNED_BYTE, data);

   BITMAPINFO bmi;
   ZeroMemory(&bmi, sizeof(BITMAPINFO));
   bmi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
   bmi.bmiHeader.biBitCount = 24;
   bmi.bmiHeader.biHeight = height;
   bmi.bmiHeader.biWidth = width;
   bmi.bmiHeader.biPlanes = 1;
   bmi.bmiHeader.biCompression = BI_RGB;
   bmi.bmiHeader.biSizeImage = size;

   Gdiplus::Bitmap bitmap(&bmi, data);
   CLSID pngClsid;
   GetEncoderClsid(L"image/png", &pngClsid);
   
   WCHAR szPath[MAX_PATH];
   GetCurrentDirectoryW(MAX_PATH, szPath);
   WCHAR szTmp[MAX_PATH];
   wsprintfW(szTmp, L"%s\\screen%05d.png", szPath, number);
   bitmap.Save(szTmp, &pngClsid);

   delete[] data;
}

// --------------------------------------------------------------------------------------------------------------------------------

#else


void screen_capture_startup()
{
   // TODO
}

void screen_capture_shutdown()
{
   // TODO
}

void screenshot(int number)
{
   // TODO
}


// --------------------------------------------------------------------------------------------------------------------------------

#endif
