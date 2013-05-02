#'
#'An R interface of taobao API
#'
#'@param method the method will be used, see detials of the method list on 
#'@param params the parameter list required in the method given above
#'@param authURL the enter point of taobao API
#'@param appKey,appSecret the key and secret of a authorized application in taobao
#'
#'
#'@return
#'it returns a \code{XMLDocument} object if connect is ok, otherwise NULL is returned.
#'
#'@seealso
#'\code{\link{getNodeSet}}
#'
#'
#'@examples
#'require(digest)
#'require(RCurl)
#'require(XML)
#'taobaoAPI("taobao.user.get",list(fields="user_id,uid,nick,sex,buyer_credit",nick="sandbox_c_4"),
#'authURL="http://gw.api.tbsandbox.com/router/rest",appKey="test",appSecret="test")
#'
#'http://s.taobao.com/search?q=e&mini=1
#'

taobaoAPI<-function(method="",params=list(),authURL,appKey,appSecret){
  apiParams=params
  apiParams$format="xml"
  apiParams$sign_method="md5"
  apiParams$app_key=appKey
  apiParams$v="2.0"
  apiParams$timestamp=as.character(Sys.time())
  
  apiParams$method=method
  strparam=createSignAndStr(apiParams,appSecret)
  url=paste(authURL,"?",strparam,sep="")
  
  txt=try(getURL(url))
  if(inherits(txt,"try-error")){
    txtree=NULL
  }else{
    txtree=xmlParse(txt)
  }
  
  return(txtree)
}

createSignAndStr<-function(apiParams,appSecret){
  
  strparam=""
  sign=appSecret
  apiParams=apiParams[order(names(apiParams))]
  np=length(apiParams)
  pnames=names(apiParams)
  for(i in 1:np){
    if(pnames[i]!="" & apiParams[[i]]!=""){
      sign=paste(sign,pnames[i],apiParams[[i]],sep="")
      strparam=paste(strparam,pnames[i],"=",URLencode(apiParams[[i]],TRUE),"&",sep="")
    }
  }
  sign=paste(sign,appSecret,sep="")
  sign=toupper(digest(sign,algo="md5",serialize=FALSE))
  strparam=paste(strparam,"sign=",sign,sep="")
  return(strparam)
}