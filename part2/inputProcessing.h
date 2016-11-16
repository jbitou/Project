#define MAX_LINE 1000

typedef struct config_info *pinfo; 
typedef struct config_info {
	int k;
	int num_of_hash;
	int L;
	int fraction;
	int iterations;
}conf_info;

int command_processing(int);
pinfo get_config_info(FILE *, int);
char *inputString(FILE *, size_t);
