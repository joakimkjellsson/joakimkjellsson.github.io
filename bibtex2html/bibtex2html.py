import bibtexparser
import bibtexparser.middlewares as bm

bibfile = "test.bib" 

# Reads a bib file and parses it
# SeparateCoAuthors means the author list will be a python list
library = bibtexparser.parse_file('test.bib', \
                                  append_middleware=[bm.SeparateCoAuthors(), bm.SplitNameParts()])

print(f"Parsed {len(library.blocks)} blocks, including:"
      f"\n\t{len(library.entries)} entries"
      f"\n\t{len(library.comments)} comments"
      f"\n\t{len(library.strings)} strings and"
      f"\n\t{len(library.preambles)} preambles")
      
def get_entry_by_author(library, lastname="Kjellsson"):
    
    entries = []
    
    for entry in library.entries:
        
        # get authors 
        authors = entry.get('author').value
        
        # get title
        title = entry.get('title').value
        
        for author in authors:
            
            # first and last names are in a list
            # i will only take the first
            last = author.last[0]
            first = author.first[0]
            
            # if lastname matches, keep it
            if last == lastname:
                entries.append(entry)
                print(" === Found %s in %s " % (last,title)) 
        
    return entries
        


def make_entries_to_html(entries):
    
    for entry in entries:
        
        authors = entry.get('author').value
        
        title = entry.get('title').value
        
        doi = entry.get('doi').value
        
        journal = entry.get('journal').value
        
        year = entry.get('year').value
        
        # initialise empty string
        html_string = ""
        
        html_string += '<p class="publication">'
        
        author_string = ""
        for author in authors:
            
            il = 0
            for last in author.last: 
                
                author_string += last
                il += 1
                
                if il == len(author.last):
                    author_string += ", "
            
            il = 0    
            for first in author.first:
                
                author_string += "%s." % (first[0],)
                il += 1
                
                if il == len(author.first):
                    author_string += ", "
                
        print(author_string)

entries = get_entry_by_author(library)

make_entries_to_html(entries)
