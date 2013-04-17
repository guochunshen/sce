\name{sendmail}
\alias{sendmail}
\title{Sends an Email with Status Information}
\usage{
  sendmail(recipient, subject = "Notification from R",
    message = "Calculation finished!")
}
\arguments{
  \item{recipient}{A valid email address}

  \item{subject}{Subject of the email message}

  \item{message}{Body of the email message}
}
\description{
  Sends an Email with Status Information
}
\details{
  This is a sample status notification way by sending email
  from \code{mail} package.
}
\note{
  To avoid spam, you have a maximum of 20 email
  notifications per day, that you can send from the same
  IP-address.
}
\author{
  Dr. Lin Himmelmann <Rmail
}

