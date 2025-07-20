# assemblmore
Automatically merges contigs provided from next-gen assemblers such as flye and canu via spanning ultra-long reads.
Future iterations will attempt to 'un-collapse' complex genomic regions such as rRNA arrays, telomeres, etc.

## NOTE
All testing was done with relatively good datasets (N50 >= 45kb) on hermaphroditic or inbred worms. 

When adding /src/ to path, make sure to NOT use ~/ shell expansion. Use the $HOME value instead.
