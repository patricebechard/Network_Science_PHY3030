# -*- coding: utf-8 -*-
"""
Crawler

PHY3030 - Projet de fin d'études
Network Science

Patrice Béchard
janvier 2017


CODE EDITED FROM : http://www.netinstructions.com/how-to-make-a-web-crawler-in-under-50-lines-of-python-code/
"""

#----------------------------------MODULES-----------------------------
from html.parser import HTMLParser  
from urllib.request import urlopen  
from urllib import parse
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import copy
from retrying import retry

print("Execution Start Time :",time.asctime())
start=time.time()            #time starts at beginning of execution
#-------------------------CLASSES-----------------------------------

class TimeError(Exception):
    """Algorithm has stalled"""
    pass

# We are going to create a class called LinkParser that inherits some
# methods from HTMLParser which is why it is passed into the definition
class LinkParser(HTMLParser):
    
    # This is a function that HTMLParser normally has
    # but we are adding some functionality to it
    def handle_starttag(self, tag, attrs):
        # We are looking for the begining of a link. Links normally look
        # like <a href="www.someurl.com"></a>
        if tag == 'base':
            for (key,value) in attrs:
                if key == 'href':
                    self.baseUrl=value
        if tag == 'a':
            for (key, value) in attrs:
                if key == 'href':
                    # We are grabbing the new URL. We are also adding the
                    # base URL to it. For example:
                    # www.netinstructions.com is the base and
                    # somepage.html is the new URL (a relative URL)
                    #
                    # We combine a relative URL with the base URL to create
                    # an absolute URL like:
                    # www.netinstructions.com/somepage.html
                    for i in range(len(value)-6):
                        if value[i:i+6]=='mailto':
                           value=''
                    newUrl = parse.urljoin(self.baseUrl, value)
                    if umontreal(newUrl)==True:
                    # And add it to our colection of links:
                        self.links = self.links + [newUrl]

    # This is a new function that we are creating to get links
    # that our spider() function will call
    def getLinks(self, url):
        self.links = []
        # Remember the base URL which will be important when creating
        # absolute URLs
        self.baseUrl=url
        # Use the urlopen function from the standard Python 3 library
        
        response=try_reading_response(url)

        # Make sure that we are looking at HTML and not other things that
        # are floating around on the internet (such as
        # JavaScript files, CSS, or .PDFs for example)
        #if response.getheader('Content-Type')=='text/html':
        try:
            htmlBytes = response.read()
            # Note that feed() handles Strings well, but not bytes
            # (A change from Python 2.x to Python 3.x)
            htmlString = htmlBytes.decode("utf-8")
            self.feed(htmlString)
            return htmlString, self.links
        #else:
        except:
            return "",[]

#------------------------------FUNCTIONS---------------------------------
def umontreal(url):
    """We limit our search to this domain only"""
    for i in range(len(url)):   
        if url[i:i+13]=='craq-astro.ca':#url[i:i+28]=='www.exoplanetes.umontreal.ca' or \
           #url[i:i+20]=='www.lps.umontreal.ca' or \
           #url[i:i+13]=='craq-astro.ca'
           #url[i:i+18]=='.phys.umontreal.ca':
            return True
    return False         
    
def delete_duplicates(links,reference,addindex=0):
    """We make sure we delete as much junk as possible from the link list"""
    
    listPointing=[]
    i=0
    temp=copy.deepcopy(links)
    for link in temp:
        i+=1
        if link in reference:
            listPointing.append(reference.index(link)+addindex)
            links.remove(link)
    
    return listPointing,links

                    
    
@retry(stop_max_delay=10)            #wait max 10 seconds
def try_reading_response(url):
    """"""
    response = urlopen(url)
    return response       
    raise TimeError("Error, time is too long")

# And finally here is our spider. It takes in an URL
# and the number of pages to search through before giving up
def crawler(start, maxPages):  
    f=open('list_of_urls.txt','w')
    fileLinks=open('fileLinks.txt','w') 
    numberNodes=open('numberNodes.txt','w') #total number of nodes (easier for plotting)
    maxNodes=1                  #node for url
    pagesToVisit = [start]        #in queue to be visited
    pagesVisited=[]             #we already visited them
    
    # The main loop. Create a LinkParser and get all the links on the page.
    # In our getLinks function we return the web page
    # (this is useful for searching for the word)
    # and we return a set of links from that web page
    # (this is useful for where to go next)
    while (len(pagesVisited)) <= maxPages and pagesToVisit != []:  #still things to visit and under the limit
        url=pagesToVisit[0]
        pagesToVisit=pagesToVisit[1:]
        print("maxNodes : ",maxNodes)
        links=[]
        pointingTo=[]         #pages where the current node points
        # Start from the beginning of our collection of pages to visit:
        
        while url in pagesVisited:          #making sure we haven't went to this page
            if pagesToVisit==[]: break
            print('-------------',url)
            print('yo')
            url=pagesToVisit[0]
            pagesToVisit=pagesToVisit[1:]
        
        f.write('{0}{1}{2}'.format(len(pagesVisited)+1,url,'\n')) #writing the list of site we've already visited
        try:
#        if True:    
            print(len(pagesVisited)+1, url)
            parser = LinkParser()
            data, links = parser.getLinks(url)  #get links from url
            links=list(set(links))              #removing duplicates 
            for link in links:
                if link == url:
                    links.remove(link)          #if link points to itself, delete from links
            newPointing,links=delete_duplicates(links,pagesVisited) 
            pointingTo+=(newPointing)
            newPointing,links=delete_duplicates(links,pagesToVisit,len(pagesVisited))
            pointingTo+=(newPointing)    
            
            #at this point, elems in links are only new links
            #add the remaining links to the list of items where the node points
            pointingTo+=([x for x in range(maxNodes,maxNodes+len(links))])
            #print(pointingTo)
            maxNodes+=len(links)            #update maxNodes

            # Add the pages that we visited to the end of our collection
            # of pages to visit:
            pagesToVisit = pagesToVisit + links
        
        except :
            #even if it failed, we write in visited because we dont want another visit
            print(" **Failed!**")  
            pointingTo=[]           #points to nothing   
        
        fileLinks.write(str(pointingTo)) 
        fileLinks.write('\n')
        pagesVisited+=[url]
        print(maxNodes,len(pagesVisited)+len(pagesToVisit))
            
    if pagesToVisit==[]:
        print('pagesToVisit is empty!')    
    
    print('Nombre de pages visitées : %d'%(len(pagesVisited)))
    print('Nombre de pages toujours en file : %d'%len(pagesToVisit))  
    numberNodes.write(str(maxNodes))
    return 

#---------------------------------MAIN----------------------------------
#crawler('http://www.exoplanetes.umontreal.ca/',10000)
#crawler('http://www.lps.umontreal.ca/fr',10000)
#crawler('http://phys.umontreal.ca/accueil/',100000)
crawler('http://craq-astro.ca/',100000)

totaltime=time.time()-start
print("Total time : %d h %d m %d s"%(totaltime//3600,totaltime//60,totaltime%60))