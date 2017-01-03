#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inputProcessing.h"

int command_processing(int argc) {
	/**Parameter -complete may be given alone at the end**/
	if (argc > 6) {
		printf("Too many arguments. Try again.\n");
		return -1;
	}
	else if (argc < 5) {
		printf("Too few arguments. Try again.\n");
		return -1;
	}
	return 0;	
}
