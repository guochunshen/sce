\name{taobaoAPI}
\alias{taobaoAPI}
\title{An R interface of taobao API}
\usage{
  taobaoAPI(method = "", params = list(), authURL, appKey,
    appSecret)
}
\arguments{
  \item{method}{the method will be used, see detials of the
  method list on}

  \item{params}{the parameter list required in the method
  given above}

  \item{authURL}{the enter point of taobao API}

  \item{appKey,appSecret}{the key and secret of a
  authorized application in taobao}
}
\value{
  it returns a \code{XMLDocument} object if connect is ok,
  otherwise NULL is returned.
}
\description{
  An R interface of taobao API
}
\examples{
require(digest)
require(RCurl)
require(XML)
taobaoAPI("taobao.user.get",list(fields="user_id,uid,nick,sex,buyer_credit",nick="sandbox_c_4"),
authURL="http://gw.api.tbsandbox.com/router/rest",appKey="test",appSecret="test")

http://s.taobao.com/search?q=e&mini=1
}
\seealso{
  \code{\link{getNodeSet}}
}

