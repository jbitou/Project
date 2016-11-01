#include <stdio.h>

int main (int argc, char **argv)
{
	FILE *fp;
	int i, v;
	char t[100];
	double d;	
	
	for (i=1; i < (argc-1); i++)
	{
		fp = fopen(argv[i],"r");
		while (fscanf(fp,"%s",t) != EOF)
		{
			if (strcmp(t,"distanceTrue:") == 0)
			{
				if ((strcmp(argv[argc-1],"-E") == 0) || (strcmp(argv[argc-1],"-C") == 0))	//Euclidean or Cosine
				{
					fscanf(fp,"%lf",&d);
					if (d > 2) 
					{
						printf("Error\n");
						return -1;
					}
				}
				else if ((strcmp(argv[argc-1],"-H") == 0) || (strcmp(argv[argc-1],"-M") == 0)) 	//Hamming or Matrix
				{
					fscanf(fp,"%d",&v);
					if (strcmp(argv[argc-1],"-H") == 0)
					{
						if (v > 23) 
						{
							printf("Error\n");
							return -1;
						}
					}
					else
					{
						if (v > 10) 
						{
							printf("Error\n");
							return -1;
						}
					}
				}
			}
		}
		fclose(fp);
		printf("Done\n");
	}
	
	return 0;
}
	
