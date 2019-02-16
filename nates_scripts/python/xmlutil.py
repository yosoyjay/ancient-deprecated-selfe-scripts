import xml.xpath as xpath
import xml
from xml.dom.ext.reader import Sax2

def ParseXML(xmltext):
  reader = Sax2.Reader()
  d = reader.fromString(xmltext)
  context = xml.xpath.Context.Context(d)
  #context.setNamespaces(namespaces)
  query = "/"
  e = xml.xpath.Compile(query)
  result = e.evaluate(context)
  node = result[0]
  return node

def ExtractValueFromNode(node):
  if node.nodeName == "#text":
    # text node
    value = node.nodeValue
  elif node.nodeValue == None:
    # element node
    value = node.localName
  else:
    # attribute node, or something we don't want
    value = node.nodeValue
  v = value.strip()
  v = v.encode('utf-8')
  return v
 
def Xpath(xpath, node):
  if type(xpath) in (str, unicode):
    xpath = xml.xpath.Compile(xpath)
  c = xml.xpath.Context.Context(node)

  nodes = xpath.evaluate(c)
  return nodes

