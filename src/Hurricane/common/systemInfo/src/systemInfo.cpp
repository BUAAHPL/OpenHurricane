/*!
 * \file systemInfo.cpp
 * \brief Main subroutines for the system info.
 * \author Rao Sihang
 * \version V2.0.0
 * \date 2022.05.02
 *
 * OpenHurricane: Open parts of Hurricane project (Highly Universal Rocket & Ramjet sImulation Codes for ANalysis and Evaluation)
 * \copyright Copyright (C) 2019-2023, Prof. Xu Xu's group at Beihang University.
 *
 * License
 *		This file is part of OpenHurricane
 *
 *		OpenHurricane is free software: you can redistribute it and/or modify it
 *		under the terms of the GNU General Public License as published by
 *		the Free Software Foundation, either version 3 of the License, or
 *		(at your option) any later version.
 *
 *		OpenHurricane is distributed in the hope that it will be useful, but WITHOUT
 *		ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *		FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *		for more details.
 *
 *		You should have received a copy of the GNU General Public License
 *		along with OpenHurricane.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 */
#include "systemInfo.hpp"
#include <cstring>
#include <iostream>

namespace OpenHurricane {

#if defined(_WIN32)

#include <Windows.h>

    void getOsInfo(std::string &osName, std::string &osVersion) {
        //先判断是否为win8.1或win10
        typedef void(__stdcall * NTPROC)(DWORD *, DWORD *, DWORD *);
        HINSTANCE hinst = LoadLibrary("ntdll.dll");
        DWORD dwMajor, dwMinor, dwBuildNumber;
        NTPROC proc = (NTPROC)GetProcAddress(hinst, "RtlGetNtVersionNumbers");
        proc(&dwMajor, &dwMinor, &dwBuildNumber);
        if (dwMajor == 6 && dwMinor == 3) //win 8.1
        {
            osName = "Microsoft Windows 8.1";
            osVersion = std::to_string(dwMajor) + ".";
            osVersion += std::to_string(dwMinor);
            return;
        }
        if (dwMajor == 10 && dwMinor == 0) //win 10
        {
            osName = "Microsoft Windows 10";
            osVersion = std::to_string(dwMajor) + ".";
            osVersion += std::to_string(dwMinor);
            return;
        }
        //判断win8.1以下的版本
        SYSTEM_INFO info;     //用SYSTEM_INFO结构判断64位AMD处理器
        GetSystemInfo(&info); //调用GetSystemInfo函数填充结构
        OSVERSIONINFOEX os;
        os.dwOSVersionInfoSize = sizeof(OSVERSIONINFOEX);
#pragma warning(disable : 4996)
        if (GetVersionEx((OSVERSIONINFO *)&os)) {
            osVersion = std::to_string(os.dwMajorVersion) + ".";
            osVersion += std::to_string(os.dwMinorVersion);

            //下面根据版本信息判断操作系统名称
            switch (os.dwMajorVersion) //判断主版本号
            {
            case 4:
                switch (os.dwMinorVersion) { //判断次版本号
                case 0:
                    if (os.dwPlatformId == VER_PLATFORM_WIN32_NT)
                        osName = "Microsoft Windows NT 4.0"; //1996年7月发布
                    else if (os.dwPlatformId == VER_PLATFORM_WIN32_WINDOWS)
                        osName = "Microsoft Windows 95";
                    break;
                case 10:
                    osName = "Microsoft Windows 98";
                    break;
                case 90:
                    osName = "Microsoft Windows Me";
                    break;
                }
                break;
            case 5:
                switch (os.dwMinorVersion) { //再比较dwMinorVersion的值
                case 0:
                    osName = "Microsoft Windows 2000"; //1999年12月发布
                    break;
                case 1:
                    osName = "Microsoft Windows XP"; //2001年8月发布
                    break;
                case 2:
                    if (os.wProductType == VER_NT_WORKSTATION &&
                        info.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_AMD64)
                        osName = "Microsoft Windows XP Professional x64 Edition";
                    else if (GetSystemMetrics(SM_SERVERR2) == 0)
                        osName = "Microsoft Windows Server 2003"; //2003年3月发布
                    else if (GetSystemMetrics(SM_SERVERR2) != 0)
                        osName = "Microsoft Windows Server 2003 R2";
                    break;
                }
                break;
            case 6:
                switch (os.dwMinorVersion) {
                case 0:
                    if (os.wProductType == VER_NT_WORKSTATION)
                        osName = "Microsoft Windows Vista";
                    else
                        osName = "Microsoft Windows Server 2008"; //服务器版本
                    break;
                case 1:
                    if (os.wProductType == VER_NT_WORKSTATION)
                        osName = "Microsoft Windows 7";
                    else
                        osName = "Microsoft Windows Server 2008 R2";
                    break;
                case 2:
                    if (os.wProductType == VER_NT_WORKSTATION)
                        osName = "Microsoft Windows 8";
                    else
                        osName = "Microsoft Windows Server 2012";
                    break;
                }
                break;
            default:
                osName = "Unknown System";
            }
        } else {
            osName = "Failed to get System info";
            osVersion = "NULL";
        }
    }

#if _MSC_VER >= 1400
#include <intrin.h>
#endif
#if defined(_WIN64)

#else
#if _MSC_VER < 1600
    void __cpuidex(INT32 CPUInfo[4], INT32 InfoType, INT32 ECXValue) {
        if (NULL == CPUInfo)
            return;
        __asm {
			mov edi, CPUInfo;
			mov eax, InfoType;
			mov ecx, ECXValue;
			cpuid;
			mov[edi], eax;
			mov[edi + 4], ebx;
			mov[edi + 8], ecx;
			mov[edi + 12], edx;
        }
    }
#endif
#if _MSC_VER < 1400
    void __cpuid(INT32 CPUInfo[4], INT32 InfoType) {
        __cpuidex(CPUInfo, InfoType, 0);
    }
#endif
#endif

    int getCpuVendor(char *pvender) {
        INT32 dwBuf[4];
        if (NULL == pvender)
            return 0;
        __cpuid(dwBuf, 0);
        *(INT32 *)&pvender[0] = dwBuf[1];
        *(INT32 *)&pvender[4] = dwBuf[3];
        *(INT32 *)&pvender[8] = dwBuf[2];
        pvender[12] = '\0';
        return 12;
    }
    int getCpuBrand(char *pbrand) {
        INT32 dwBuf[4];
        if (NULL == pbrand)
            return 0;
        __cpuid(dwBuf, 0x80000000);
        if (dwBuf[0] < 0x80000004)
            return 0;
        __cpuid((INT32 *)&pbrand[0], 0x80000002);
        __cpuid((INT32 *)&pbrand[16], 0x80000003);
        __cpuid((INT32 *)&pbrand[32], 0x80000004);
        pbrand[48] = '\0';
        return 48;
    }

    void getCpuInfo(int &np, std::string &cpuVendor, std::string &cpuBrand) {
        char szBuf[64];

        SYSTEM_INFO si;
        GetSystemInfo(&si);
        np = si.dwNumberOfProcessors;
        //cout << si.dwNumberOfProcessors << endl;
        getCpuVendor(szBuf);
        cpuVendor = szBuf;
        //cout << szBuf << endl;
        getCpuBrand(szBuf);
        cpuBrand = szBuf;
        //cout << szBuf << endl;
    }

#endif //_WIN32

#if defined(__linux__) || defined(LINUX)
#include <sys/sysinfo.h>
#include <sys/utsname.h>
    void getOsInfo(std::string &osName, std::string &osVersion) {
        struct utsname u;
        uname(&u);
        osName = u.sysname;
        osVersion = u.version;
    }

    char *my_strstr(const char *s1, const char *s2) {
        const char *p = s1;
        const size_t len = strlen(s2);
        int i = 0;
        for (; (p = strchr(p, *s2)) != 0; p++) {
            if (strncmp(p, s2, len) == 0) {
                return (char *)p;
            }
        }
        return (0);
    }

    void getCpuInfo(int &np, std::string &cpuVendor, std::string &cpuBrand) {
        np = get_nprocs();
        FILE *fp = fopen("/proc/cpuinfo", "r");
        if (NULL == fp)
            printf("failed to open cpuinfo\n");
        char szTest[1000] = {0};
        // read file line by line
        char *p;
        while (!feof(fp)) {
            memset(szTest, 0, sizeof(szTest));
            if (fgets(szTest, sizeof(szTest) - 1, fp) != nullptr) {
                // leave out \n
            }
            p = my_strstr(szTest, "vendor_id   :");
            if (p != NULL)
                cpuVendor = p;
            p = my_strstr(szTest, "model name  :");
            if (p != NULL)
                cpuBrand = p;
            //printf("%s", szTest);
        }
        fclose(fp);
    }
#endif // __linux__

} // namespace OpenHurricane
