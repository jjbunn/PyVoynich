# -*- coding: utf-8 -*-
'''
Created on May 26, 2015

Reads a file of text, e.g. the EVA or Voyn_101 transcriptions of the Voynich manuscript, or some other text file and
generates various statistics from it.

Example usage in main. 

@author: Julian Bunn, jjbunn@gmail.com

'''
import math
from collections import Counter

def myrange(a,b):
    return list(range(a,b))


V101_Gallows = 'fghjklruvFGHJKLRUV'
EVA_Gallows = 'fkptFKPT'

# from http://www.voynich.nu/layout.html
Quires = {'Q1': myrange(1,9), 'Q2': myrange(9,17), 'Q3': myrange(17,25), 'Q4': myrange(25,33), 'Q5': myrange(33,41),
          'Q6': myrange(41,49), 'Q7': myrange(49,57), 'Q8': myrange(57,67), 'Q9': myrange(67,69), 'Q10': myrange(69,71),
          'Q11': myrange(71,73), 'Q12': myrange(73,75), 'Q13': myrange(75,85), 'Q14': myrange(85,87),
          'Q15': myrange(87,91), 'Q16': myrange(91,93), 'Q17': myrange(93,97), 'Q18': myrange(97,99),
          'Q19': myrange(99,103), 'Q20': myrange(103,117)}

Sections = {'Herbal': myrange(1,57) + [87, 90, 93, 94, 95, 96],
            'Astrological': myrange(67,75),
            'Balneological': myrange(75,85),
            'Cosmological': [85,86],
            'Pharmaceutical': [88,89,99,100,101,102],
            'Recipes': myrange(103,117),
            'Mixed': myrange(57,67)}

LanguageHands = {'A1': myrange(1,26) + myrange(27,31) + [32,35,36,37,38,42,44,45,47,49,51,52,53,54,56],
                 'B2': [26,31,33,34,39,40,41,43,46,48,50,55,57] + myrange(75,85),
                 'A4': [87,88,93,96,99,100,101],
                 'BX': [103,104,105,106],
                 'B': [66] + myrange(107,117),
                 'A': [58,89,90,102],
                 'B?': myrange(67,75)}

Special = {'HerbalRecipeAB': myrange(107,117) + myrange(1,26),
           'HerbalAB': myrange(1,57) + [87,90] + myrange(93,97),
           'HerbalBalneoAB': myrange(1,26) + myrange(75,85),
           'HerbalFakeAB': myrange(13,26) + myrange(1,13),
           'HerbalAstroAB': myrange(1,13) + myrange(67,75),
           'PharmaRecipeAB': [88,89,99,100,101,102] + myrange(103,117),
           'AllAB': myrange(1,117),
           'HerbalA': myrange(1,26),
           'HerbalB': [26,31,33,34,39,40,41,43,46,48,50,55],
           'JustA': myrange(1,26) + myrange(27,31) + [32,35,36,37,38,42,44,45,47,49,51,52,53,54,56] + [87,88,93,96,99,100,101] + [58,89,90,102],
           'JustB': [26,31,33,34,39,40,41,43,46,48,50,55,57] + myrange(75,85) + [103,104,105,106] + [66] + myrange(107,117) + myrange(67,75),
           }

basechars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789?!*'

ZodiacNames = {'70v2': 'Pisces',
               '70v1': 'Aries (Dark)',
               '71r': 'Aries (Light)',
               '71v': 'Taurus (Light)',
               '72r1': 'Taurus (Dark)',
               '72r2': 'Gemini',
               '72r3': 'Cancer',
               '72v3': 'Leo',
               '72v2': 'Virgo',
               '72v1': 'Libra',
               '73r': 'Scorpio',
               '73v': 'Sagittarius',
               '67r1': '67r1',
               '67r2': '67r2',
               '68r1': '68r1',
               '68r2': '68r2',
               '68r3': '68r3',
               '68v2': '68v2',
               
               }
StarCounts = { # my counting ...
'103r':19,
'103v':14,
'104r':13,
'104v':13,
'105r':10,
'105v':10,
'106r':16,
'106v':14,
'107r':15,
'107v':15,
'108r':16,
'108v':16,
# missing 109r,109v,110r,110v
'111r':17,
'111v':19,
'112r':12,
'112v':13,
'113r':16,
'113v':15,
'114r':13,
'114v':12,
'115r':13,
'115v':13,
'116r':10
}


class Statistics(object):

    def __init__(self, inputfile, transcription='EVA', voynich=True, lineseperator='|'):
        print ('Processing file', inputfile)
        linecount = 0
        nfolio = 1

        split_types = ['diagram', 'figure', 'plant', 'pond', 'star', ]
        
        self.basechars = basechars
        if transcription == 'EVA':
            self.Gallows = EVA_Gallows
        else:
            self.Gallows = V101_Gallows
            
        self.Quires = Quires
        self.Sections = Sections
        self.LanguageHands = LanguageHands
        self.Special = Special
        
        self.foliolines = {}
        self.foliofeatures = {}
        self.folioparagraphs = {}
        self.folioendlines = {}
        self.foliolabels = {}
        self.foliostars = {}
        self.all_text = ''
        
        last_folio = ''
        last_paragraph = ''
        line_num_on_folio = 0
        paragraphs = []

        try:
            for line in open(inputfile, encoding="ISO-8859-1").readlines():
                #print (line)
                linecount += 1
                tokens = line.split('>')
                folioline = tokens[0].strip('<').split('.')
                folio = folioline[0][1:]

                if folio != last_folio:
                    if last_folio != '':
                        self.folioparagraphs[last_folio] = paragraphs
                    line_num_on_folio = 0
                    last_folio = folio
                    last_paragraph = ''
                    paragraphs = []
                line_num_on_folio += 1

                if voynich:
                    text = tokens[1]
                else:
                    text = tokens[0]

                features = []

                if voynich:
                    if transcription == 'EVA':
                        linetype = folioline[1]
                 
                        line_start_indices = [0]
                        text = text.strip()
                        text = text.replace('!','.')

                        if folioline[1][0] == 'P':
                            if last_paragraph != folioline[1]:
                                paragraphs.append(line_num_on_folio)
                            last_paragraph = folioline[1]
                        if "{" in text:
                            # remove curly brace keywords  

                            while "{" in text:
                                feature = text[text.index('{')+1:text.index('}')]

                                position = text.index('{')
                                if position > 0 and feature in split_types:
                                    features.append( (feature, position) )
                                text = text[:position] + text[text.index('}')+1:]
                                        

                        endparagraph = text.endswith('=')
                        text = text.replace('.',' ')
                        text = text.replace('-',' ')
                        text = text.replace('=','')

                        if linetype[0] == 'S' or \
                            (folio == '67r2' and linetype[0] in 'LXZ') or \
                            (folio == '67v2' and linetype[0] == 'F') or \
                            (folio == '67v1' and linetype[0] == 'X') or \
                            (folio == '68r3' and linetype[0] == 'X') or \
                            (folio == '68v1' and linetype[0] == 'X'):
                            
                            star = [text]
                            
                            if folio == '68v2':
                                star = text.split()
                            for s in star:
                                if folio in self.foliostars:
                                    self.foliostars[folio].append(s)
                                else:
                                    self.foliostars[folio] = [s]
                     
                        labels = False
                    else:
                        linenumber = folioline[1]
                        labels = not self.isInteger(linenumber)
                        text = self.filterV101(text)
                        endparagraph = text.endswith('=')
                        text = text.replace('.',' ')
                        text = text.replace('-','')
                        text = text.replace('=','')
                else: # Voyn_101
                    text = self.filterNatural(line)
                    if len(text) == 0:
                        continue
                    folio = nfolio
                    endparagraph = False
                    labels = False
                
                self.all_text += text + ' '    
                if folio in self.foliolines:
                    self.foliolines[folio].append((text)) 
                    self.foliofeatures[folio].append((features)) 

                else:
                    self.foliolines[folio] = [text] 
                    self.foliofeatures[folio] = [features] 
               
                if endparagraph:
                    if folio in self.folioendlines:
                        self.folioendlines[folio] += lineseperator + text 
                    else:
                        self.folioendlines[folio] = '' + text 
                if labels:
                    if folio in self.foliolabels:
                        self.foliolabels[folio] += lineseperator + text
                    else:
                        self.foliolabels[folio] = '' + text  
             
            self.folioparagraphs[folio] = paragraphs          
            #self.foliostars['70v2'].append('otylal') # star in centre of chart, next to upper fish
            #self.foliostars['72r2'].append(' ') # unlabelled star 
        except Exception as err:
            print ('Error', err)

        print ('Total of', linecount, 'lines in', inputfile) 
           
    def gallowsWord(self, word):
        for char in self.Gallows:
            if word.count(char) > 0: return True
        return False
    
    def getSection(self, folio):
        number = folio
        if folio == 'rose': return 'rose'
        while (not number.isdigit() and len(number) > 0):
            number = number[:len(number)-1]
        for section in self.Sections.keys():
            if int(number) in self.Sections.get(section): return section
        return 'Unknown'
    
    def isSpecial(self,key,folio):
        number = folio
        if folio == 'rose': return False
        while (not number.isdigit() and len(number) > 0):
            number = number[:len(number)-1]
        abnumbers = self.Special.get(key)
        return int(number) in abnumbers
        
    def getLanguageHand(self,folio):
        number = str(folio)
        if folio == 'rose': return 'rose'
        while (not number.isdigit() and len(number) > 0):
            number = number[:len(number)-1]
        for languagehand in self.LanguageHands.keys():
            if int(number) in self.LanguageHands.get(languagehand): return languagehand
        return 'Unknown'
    
    def getQuire(self,folio):
        number = folio
        if folio == 'rose': return 'rose'
        while (not number.isdigit() and len(number) > 0):
            number = number[:len(number)-1]
        for quire in self.Quires.keys():
            if int(number) in self.Quires.get(quire): return quire
        return 'Unknown'
    
    def filterNatural(self,text):
        # regularises natural language texts in e.g. English
        text = text.lower()
        text = text.strip(' ')
        text = text.replace('\n',' ')
        text = text.replace('\r',' ')
        text = text.replace(',',' ')
        text = text.replace('.',' ')
        text = text.replace(':',' ')
        text = ' '.join(text.split()) # remove multiple spaces
        return text
           
    def filterV101(self,line):
        # takes a line of Voyn101 encoding and replaces some rarer glyphs with their likely more common values
        # for very rare glyphs, replaces them by "?"
        line = line.strip(' \r\n')
        line = line.replace(',','')
        
        replacements = {'.':' ', ',':'', '3':'2',  '5':'2',  '+':'2',  '%':'2',  '#':'2',  '6':'8',  '7':'8',  'A':'a',  
                        'X':'y',  'I':'ii',  'C':'cc',  'z':'iy',  'Z':'iiy',  'j':'g',  'u':'g',  
                        'd':'ccc',  'U':'G',  'P':'ip',  'Y':'y',  '$':'s',  'S':'cs',  't':'s',  
                        'q':'iip',  'm':'iiN',  'M':'iiiN',  'n':'iN',  'Y':'y',  '!':'2',  ')':'9',  
                        '*':'p',  'b':'y',  'J':'G',  'E':'c',  'x':'y',  'B':'cccN',  'D':'ccN',  
                        'T':'s',  'Q':'p',  'W':'g',  'w':'g',  'V':'G',  '&':'8',  
                        #'f':'g', 'h':'g',  'k':'g',  'l':'r',  'F':'G',  'H':'G',  'K':'G',  'L':'r',  'R':'r',
                        }
         
        for seek, repl in replacements.items():
            line = line.replace(seek,repl)
        
        rare = set(line) - set(self.basechars+' ')    
        for c in rare:
            line = line.replace(c,'?')
                
        return line
            
    def isInteger(self,s):
        try:
            int(s)
            return True
        except ValueError:
            return False
    
    def getNgrams(self, string, MINGRAM, NGRAM, ngrams):
        # adds all ngrams of length MINGRAM to NGRAM in string to dict ngrams
        for i in range(0,len(string)):
            for n in range(MINGRAM-1,NGRAM):
                if i+n < len(string):
                    ngram = string[i:i+n+1]
                    if ngram.count('_') == 0: continue
                    if ngram.endswith('_') or ngram.startswith('_'): continue
                    if ngram in ngrams:
                        ngrams[ngram] = ngrams.get(ngram) + 1 
                    else:
                        ngrams[ngram] = 1 
        return ngrams

    def entropy(self,s):
        # Note that we use log base 2 for entropy in bits
        p, lns = Counter(s), float(len(s))
        return -sum( count/lns * math.log(count/lns, 2) for count in p.values())
    
    def character_counts(self,s):
        return Counter(s).items()

if __name__ == '__main__':   
    stats = Statistics('voyn_101.txt', transcription='V101', voynich=True, lineseperator='|')

    print (stats.all_text[:80])    
    print ('Voynich V101 Entropy',round(stats.entropy(stats.all_text),3))
    print ('')

    stats = Statistics('eva_takeshi.txt', voynich=True, lineseperator='|')

    print (stats.all_text[:80])    
    print ('Voynich EVA Entropy',round(stats.entropy(stats.all_text),3))
    print ('')
    
    chars = stats.basechars
    string = ''
    import random
    while len(string) < 5000:
        string += random.choice(chars) 
    print (string[:80]) 
    print ('Random Entropy',round(stats.entropy(string),3)) 
    print ('')
    
       
