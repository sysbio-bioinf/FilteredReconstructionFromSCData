/*	inserts only uniq values
*/

#ifndef _ENTRIES_H_
#define _ENTRIES_H_

#include <stdio.h>

typedef struct dep_s
{	char	color;
	const char	*s;
	int	count;
	struct dep_s	*left,*right,*parent;
} dep_t;

dep_t *dep_first(dep_t *p);
dep_t *dep_next(dep_t *p);

dep_t *dep_add(dep_t *p,const char *s);
dep_t *dep_find(dep_t *p,const char *s);

size_t dep_count(dep_t *p);
size_t dep_count_th(dep_t *p,int th);

void dep_free(dep_t *p);


typedef struct entry_s
{	char    color;
	const char	*s;
	struct dep_s	*depends;
	struct entry_s	*left,*right,*parent;
} entry_t;

entry_t *entry_first(entry_t *p);
entry_t *entry_next(entry_t *p);

entry_t *entry_add(entry_t *p,const char *s,const char *d);
entry_t *entry_find(entry_t *p,const char *s);

size_t entry_count(entry_t *p);

void entry_free(entry_t *p);

#endif
