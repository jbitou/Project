
typedef struct rating_info *ratingp;
typedef struct rating_info {
	int itemId;
	int rate;	
}rating;

typedef struct user_info {
	int userId;
	ratingp ratings;	
}user;


int command_processing(char **, int, int *, int *, int *);
int count_data(FILE *, int *, int *);
user *create_users(FILE *, int, int);
