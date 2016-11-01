
typedef struct nnr_node *nnrp;
typedef struct nnr_node
{
	char key[15];
	nnrp next;
}nnr;

void insert_nnrlist(char *, nnrp *);
void destroy_nnrlist(nnrp *, FILE *);
void print_nnrlist(nnrp);
