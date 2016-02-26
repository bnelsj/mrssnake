/*

This program was created at:  Fri Feb 26 14:47:00 2016
This program was created by:  Brad Nelson


Contact: bnelsj@gmail.com

Organization: University of Washington

The MIT License (MIT)

Copyright (c) <2016> <Brad Nelson>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


*/

#include <string>
#include <iostream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <unistd.h>



struct options{
   std::string file;
}globalOpts;

static const char *optString = "b:p:n:c";

//-------------------------------   OPTIONS   --------------------------------
int parseOpts(int argc, char** argv)
    {
    int opt = 0;
    globalOpts.file = "NA";
    opt = getopt(argc, argv, optString);
    while(opt != -1){
	switch(opt){
		case 'b':
		{
		 break;
		}
		case 'p':
		{
		 break;
		}
		case 'n':
		{
		 break;
		}
		case 'c':
		{
		 break;
		}
		case '?':
		{
		 break;
		}
	}
 
  opt = getopt( argc, argv, optString ); 
   }
return 1;
}
//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  :

 Function does   :

 Function returns:

*/
void sub()
{
}

//-------------------------------    MAIN     --------------------------------
/*
 Comments:
*/

int main( int argc, char** argv)
{
int parse = parseOpts(argc, argv);

return 0;
}
