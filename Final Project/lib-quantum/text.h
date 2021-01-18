/*
text.h
Written by Walter O. Krawec
Copyright (c) 2017

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef _WOK_TEXT_H_
#define _WOK_TEXT_H_

#include <string>
#include <vector>

namespace wok
{
    class Text
    {
    public:
        struct Pair
        {
            Pair() : start(0), end(0) {};
            Pair(char _start, char _end) : start(_start), end(_end) {};
            char start, end;
        };
        
        Text()
        {
            buffer = NULL;
            bufferSize = 0;
            textSize = 0;
            
            lastLookup = 0;
        }
        ~Text()
        {
            cleanup();
        }
        
        void reset()
        {
            ignore.clear();
            delimiter.clear();
            collection.clear();
            
            calculateTextSize();
        }
        void cleanup()
        {
            if(buffer)
                delete[] buffer;
            buffer = NULL;
            bufferSize = 0;
            textSize = 0;
        }
        
        void parse(const std::string& textToParse)
        {
            text = textToParse;
            
            if(bufferSize < text.size()+2)
            {
                if(buffer)
                    delete[] buffer;
                
                buffer = new char[text.size()+2];
                bufferSize = text.size()+2;
            }
            
            calculateTextSize();
        }
        
        /*void combineToSingle(char c)
         {
         char* temp = new char[text.size()+2];
         temp[0] = 0;
         for(unsigned int i=0; i<text.size(); ++i)
         {
         temp[i] = 0;
         }
         }*/
        
        void addIgnore(char c)
        {
            ignore.push_back(c);
            calculateTextSize();
        }
        void addDelimiter(char c)
        {
            delimiter.push_back(c);
            calculateTextSize();
        }
        
        void addCollection(char c)
        {
            collection.push_back(c);
            calculateTextSize();
        }
        
        void addCollection(Pair c)
        {
            pairCollection.push_back(c);
            calculateTextSize();
        }
        
        const std::string& operator[](unsigned int index)
        {
            if(lastLookup == index)
                return localString;
            
            get(index, localString);
            lastLookup = index;
            return localString;
        }
        void get(unsigned int index, std::string& output)
        {
            output = "";
            if(text.size() == 0 || !buffer)
                return;
            
            unsigned int counter = 0, bufferIndex=0, stringIndex=0;
            
            buffer[0] = 0;
            
            int pairIndex=0;
            
            while(counter <= index)
            {
                if(stringIndex >= text.size())
                    break;
                while(stringIndex < text.size() && !checkVector(text[stringIndex], delimiter))
                {
                    if(counter == index)
                    {
                        if(checkVector(text[stringIndex], collection))
                        {
                            ++stringIndex;
                            while(stringIndex < text.size() && bufferIndex < bufferSize-1 && !checkVector(text[stringIndex], collection))// && !checkPairVector(text[stringIndex], pairCollection, true, true, pairIndex))
                            {
                                buffer[bufferIndex] = text[stringIndex];
                                ++bufferIndex;
                                buffer[bufferIndex] = 0;
                                
                                ++stringIndex;
                            }
                        }
                        else if(checkPairVector(text[stringIndex], pairCollection, true, false, pairIndex))
                        {
                            // skip ahead to closing pair
                            ++stringIndex;
                            int depth=0;
                            while(stringIndex < text.size())
                            {
                                if(text[stringIndex] == pairCollection[pairIndex].end)
                                {
                                    if(depth == 0)
                                    {
                                        break;
                                    }
                                    --depth;
                                }
                                else if(text[stringIndex] == pairCollection[pairIndex].start)
                                    ++ depth;
                                
                                buffer[bufferIndex] = text[stringIndex];
                                ++bufferIndex;
                                buffer[bufferIndex] = 0;
                                
                                ++stringIndex;
                            }
                        }
                        else if(!checkVector(text[stringIndex], ignore) && bufferIndex < bufferSize-1)
                        {
                            buffer[bufferIndex] = text[stringIndex];
                            ++bufferIndex;
                            buffer[bufferIndex] = 0;
                        }
                        
                    }
                    else if(checkVector(text[stringIndex], collection))
                    {
                        ++stringIndex;
                        while(stringIndex < text.size() && !checkVector(text[stringIndex], collection))
                            ++stringIndex;
                    }
                    else if(checkPairVector(text[stringIndex], pairCollection, true, false, pairIndex))
                    {
                        // skip ahead to closing pair
                        ++stringIndex;
                        int depth=0;
                        while(stringIndex < text.size())
                        {
                            if(text[stringIndex] == pairCollection[pairIndex].end)
                            {
                                if(depth == 0)
                                    break;
                                --depth;
                            }
                            else if(text[stringIndex] == pairCollection[pairIndex].start)
                                ++ depth;
                            ++stringIndex;
                        }
                    }
                    ++stringIndex;
                }
                ++counter;
                ++stringIndex;
            }
            
            output = std::string(buffer);
        }
        
        unsigned int size()
        {
            return textSize;
        }
    private:
        void calculateTextSize()
        {
            textSize = 1;
            int pairIndex=0;
            for(unsigned int i=0; i<text.size(); ++i)
            {
                if(checkVector(text[i], delimiter))
                    ++textSize;
                else if(checkVector(text[i], collection))
                {
                    ++i;
                    while(i < text.size() && !checkVector(text[i], collection))
                        ++i;
                }
                else if(checkPairVector(text[i], pairCollection, true, false, pairIndex))
                {
                    ++i;
                    int depth = 0;
                    while(i < text.size())
                    {
                        if(text[i] == pairCollection[pairIndex].end)
                        {
                            if(depth == 0)
                            {
                                break;
                            }
                            --depth;
                        }
                        else if(text[i] == pairCollection[pairIndex].start)
                            ++ depth;
                        ++i;
                    }
                }
            }
            
            lastLookup = textSize+1;
        }
        unsigned int skipIgnore(unsigned int index)
        {
            bool found;
            while(index < text.size())
            {
                found = false;
                for(unsigned int i=0; i<ignore.size(); ++i)
                {
                    if(text[index] == ignore[i])
                    {
                        found = true;
                        break;
                    }
                }
                
                if(!found)
                    return index;
                
                ++index;
            }
            return index;
        }
        
        bool checkVector(char c, std::vector <char>& v)
        {
            for(unsigned int i=0; i<v.size(); ++i)
            {
                if(c == v[i])
                    return true;
            }
            return false;
        }
        bool checkPairVector(char c, std::vector <Pair>& v, bool checkStart, bool checkEnd, int& pairIndex)
        {
            for(unsigned int i=0; i<v.size(); ++i)
            {
                if((checkStart && c == v[i].start) || (checkEnd && c == v[i].end))
                {
                    pairIndex = i;
                    return true;
                }
            }
            return false;
        }
        
        std::string text;
        std::vector <char> ignore, delimiter, collection;
        std::vector <Pair> pairCollection;
        char* buffer;
        unsigned int bufferSize;
        unsigned int textSize;
        
        std::string localString;
        unsigned int lastLookup;
    };
}

#endif
