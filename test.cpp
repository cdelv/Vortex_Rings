#include <iostream>
#include<string>
#include <sstream>
#include <regex>

//#include <sys/types.h>
#include <sys/stat.h>
//#include <unistd.h>
//#include <stdio.h>
//#include <stdlib.h>

/*#ifndef _WIN32
#include <sys/stat.h>  // mkdir
#else
#include <direct.h>    // _mkdir
#define mkdir(dir, mode) _mkdir(dir)
#endif*/

int main(int argc, char const *argv[])
{
    int a = mkdir("test",0777);

    return 0;
}