/*	inserts only uniq values
*/

#include "entries.h"

#include <stdlib.h>
#include <string.h>

static const char RED=0;
static const char BLACK=1;

/************************** dep **********************/
int dep_is_red(dep_t *p)
{	if (!p)
	{	return 0;
	}
	return RED==p->color;
}

dep_t *dep_most_left(dep_t *p)
{	while (p && p->left)
	{	p=p->left;
	}
	return p;
}

dep_t *dep_up(dep_t *p,const char *s)
{	while (p && strcmp(s,p->s)>=0)
	{	p=p->parent;
	}
	return p;
}

dep_t *dep_next(dep_t *p)
{	if (p->right)
	{	return dep_most_left(p->right);
	}
	return dep_up(p,p->s);
}

dep_t *dep_first(dep_t *p)
{	return dep_most_left(p);
}

dep_t *dep_insert(dep_t *p,dep_t *n)
{	if (NULL==p)
	{	p=n;
	}
	if (strcmp(p->s,n->s)>0)
	{	p->left=dep_insert(p->left,n);
		p->left->parent=p;
	}
	else if (strcmp(p->s,n->s)<0)
	{	p->right=dep_insert(p->right,n);
		p->right->parent=p;
	}
	return p;
}

dep_t *dep_find(dep_t *p,const char *s)
{	if (!p)
	{	return p;
	}
	int discr=strcmp(p->s,s);
	if (discr==0)
	{	return p;
	}
	else if (discr>0)
	{	return dep_find(p->left,s);
	}
	else
	{	return dep_find(p->right,s);
	}
}	

dep_t *dep_sibling(dep_t *p)
{	if (p->parent)
	{	if (p==p->parent->left)
		{	return p->parent->right;
		}
		return p->parent->left;
	}
	return NULL;
}

dep_t *dep_uncle(dep_t *p)
{	if (p && p->parent && p->parent->parent)
	{	return dep_sibling(p->parent);
	}
	return NULL;
}

void dep_rotate_left(dep_t *p)
{	if (p->parent)
	{	if (p->parent->left==p)
		{	p->parent->left=p->right;
		}
		else
		{	p->parent->right=p->right;
		}
	}
	dep_t *n=p->right;
	p->right=n->left;
	n->left=p;
	n->parent=p->parent;
	p->parent=n;
	if (p->right)
	{	p->right->parent=p;
	}
}

void dep_rotate_right(dep_t *p)
{	if (p->parent)
	{	if (p->parent->left==p)
		{	p->parent->left=p->left;
		}
		else
		{	p->parent->right=p->left;
		}
	}
	dep_t *n=p->left;
	p->left=n->right;
	n->right=p;
	n->parent=p->parent;
	p->parent=n;
	if (p->left)
	{	p->left->parent=p;
	}
}

void dep_rebalance(dep_t *n)
{	if (n->parent==NULL)
	{	n->color=BLACK;
	}
	else if (!dep_is_red(n->parent))
	{
	}
	else if (dep_is_red(dep_uncle(n)))
	{	n->parent->color=BLACK;
		dep_uncle(n)->color=BLACK;
		n->parent->parent->color=RED;
		dep_rebalance(n->parent->parent);
	}
	else
	{	dep_t*  p=n->parent;
		dep_t*  g=p->parent;
		if (g->left && n==g->left->right)
		{	dep_rotate_left(p);
			n=n->left;
		}
		else if (g->right && n==g->right->left)
		{	dep_rotate_right(p);
			n=n->right;
		}

		p=n->parent;
		g=p->parent;
		if (n==p->left)
		{       dep_rotate_right(g);
		}
		else
		{	dep_rotate_left(g);
		}
		p->color=BLACK;
		g->color=RED;
	}
}

dep_t *dep_add(dep_t *p,const char *s)
{	
	dep_t *n=dep_find(p,s);
	if (n)
	{	n->count++;
		return p;
	}
	else
	{	n=malloc(sizeof(dep_t));
		n->color=RED;
		n->s=s;
		n->count=1;
		n->left=NULL;
		n->right=NULL;
		n->parent=NULL;
		dep_insert(p,n);
		dep_rebalance(n);
		while (n->parent)
		{	n=n->parent;
		}
		return n;
	}
}

void dep_free(dep_t *p)
{	if (!p)
	{	return;
	}
	if (p->left)
	{	dep_free(p->left);
	}
	if (p->right)
	{	dep_free(p->right);
	}
	free(p);
}

size_t dep_count(dep_t *p)
{	if (!p)
	{	return 0;
	}
	return dep_count(p->left)+1+dep_count(p->right);
}

size_t dep_count_th(dep_t *p,int th)
{       if (!p)
	{       return 0;
	}
	return dep_count(p->left)+(p->count>=th)?1:0+dep_count(p->right);
}

/************************** entry **********************/
int entry_is_red(entry_t *p)
{	if (!p)
	{	return 0;
	}
	return RED==p->color;
}

entry_t *entry_most_left(entry_t *p)
{	while (p && p->left)
	{	p=p->left;
	}
	return p;
}

entry_t *entry_up(entry_t *p,const char *s)
{	while (p && strcmp(s,p->s)>=0)
	{	p=p->parent;
	}
	return p;
}

entry_t *entry_next(entry_t *p)
{	if (p->right)
	{	return entry_most_left(p->right);
	}
	return entry_up(p,p->s);
}

entry_t *entry_first(entry_t *p)
{	return entry_most_left(p);
}

entry_t *entry_insert(entry_t *p,entry_t *n)
{	if (NULL==p)
	{	p=n;
	}
	if (strcmp(p->s,n->s)>0)
	{	p->left=entry_insert(p->left,n);
		p->left->parent=p;
	}
	else if (strcmp(p->s,n->s)<0)
	{	p->right=entry_insert(p->right,n);
		p->right->parent=p;
	}
	return p;
}

entry_t *entry_find(entry_t *p,const char *s)
{	if (!p)
	{	return p;
	}
	int discr=strcmp(p->s,s);
	if (discr==0)
	{	return p;
	}
	else if (discr>0)
	{	return entry_find(p->left,s);
	}
	else
	{	return entry_find(p->right,s);
	}
}	

entry_t *entry_sibling(entry_t *p)
{	if (p->parent)
	{	if (p==p->parent->left)
		{	return p->parent->right;
		}
		return p->parent->left;
	}
	return NULL;
}

entry_t *entry_uncle(entry_t *p)
{	if (p && p->parent && p->parent->parent)
	{	return entry_sibling(p->parent);
	}
	return NULL;
}

void entry_rotate_left(entry_t *p)
{	if (p->parent)
	{	if (p->parent->left==p)
		{	p->parent->left=p->right;
		}
		else
		{	p->parent->right=p->right;
		}
	}
	entry_t *n=p->right;
	p->right=n->left;
	n->left=p;
	n->parent=p->parent;
	p->parent=n;
	if (p->right)
	{	p->right->parent=p;
	}
}

void entry_rotate_right(entry_t *p)
{	if (p->parent)
	{	if (p->parent->left==p)
		{	p->parent->left=p->left;
		}
		else
		{	p->parent->right=p->left;
		}
	}
	entry_t *n=p->left;
	p->left=n->right;
	n->right=p;
	n->parent=p->parent;
	p->parent=n;
	if (p->left)
	{	p->left->parent=p;
	}
}

void entry_rebalance(entry_t *n)
{	if (n->parent==NULL)
	{	n->color=BLACK;
	}
	else if (!entry_is_red(n->parent))
	{
	}
	else if (entry_is_red(entry_uncle(n)))
	{	n->parent->color=BLACK;
		entry_uncle(n)->color=BLACK;
		n->parent->parent->color=RED;
		entry_rebalance(n->parent->parent);
	}
	else
	{	entry_t*  p=n->parent;
		entry_t*  g=p->parent;
		if (g->left && n==g->left->right)
		{	entry_rotate_left(p);
			n=n->left;
		}
		else if (g->right && n==g->right->left)
		{	entry_rotate_right(p);
			n=n->right;
		}

		p=n->parent;
		g=p->parent;
		if (n==p->left)
		{       entry_rotate_right(g);
		}
		else
		{	entry_rotate_left(g);
		}
		p->color=BLACK;
		g->color=RED;
	}
}

entry_t *entry_add(entry_t *p,const char *s,const char *d)
{	
	entry_t *n=entry_find(p,s);
	if (n)
	{	n->depends=dep_add(n->depends,d);
		return p;
	}
	else
	{	n=malloc(sizeof(entry_t));
		n->color=RED;
		n->s=s;
		n->depends=dep_add(NULL,d);
		n->left=NULL;
		n->right=NULL;
		n->parent=NULL;
		entry_insert(p,n);
		entry_rebalance(n);
		while (n->parent)
		{	n=n->parent;
		}
		return n;
	}
}

void entry_free(entry_t *p)
{	if (!p)
	{	return;
	}
	if (p->left)
	{	entry_free(p->left);
	}
	if (p->right)
	{	entry_free(p->right);
	}
	dep_free(p->depends);
	free(p);
}

size_t entry_count(entry_t *p)
{	if (!p)
	{	return 0;
	}
	return entry_count(p->left)+1+entry_count(p->right);
}
