typedef struct _Stack Stack;

Stack * stack_new (int size);

void   stack_push    (Stack * s, void * p);
void * stack_pop     (Stack * s);
void * stack_index   (Stack * s, int i);
void   stack_destroy (Stack * s);
