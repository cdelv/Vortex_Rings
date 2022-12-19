#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int main()
{   
    double Z0 = 0.0;
    double a = 0.2;
    double x = 1.0;
    double y = 1.2;
    double z = 0.9;

    char command[200];
    sprintf (command, "./integral %g %g %g %g %g", Z0, a, x, y, z);
    fprintf (stderr, " %s \n", command);

    FILE *cmd=popen(command, "r");

    char buf[200]={0x0};
    fgets(buf, sizeof(buf), cmd);
    fprintf (stderr, " %s \n", buf);

    char *token = strtok(buf, ",");
    double vals[3];

    for (int i = 0; i < 3; ++i)
    {   
        vals[i] = atof(token);
        printf("%g\n", vals[i]);
        token = strtok(NULL, ",");
    }


    pclose(cmd);



    return 0;
}