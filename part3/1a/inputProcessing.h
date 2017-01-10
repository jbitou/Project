#define STARTSIZE 10 
#define f 10
#define L 5
#define k 4

typedef struct rating_info *ratingp;
typedef struct rating_info {
	int itemId;
	double rate;	
}rating;

typedef struct user_info {
	int userId;
	ratingp ratings;	
}user;


int command_processing(char **, int, int *, int *, int *);
int *create_items(FILE *, int *, int *, int *);
user *create_users(FILE *, int *, int, int);
void quickSort(int *, int, int);
void swap(int *, int *);
int partition(int *, int, int);
