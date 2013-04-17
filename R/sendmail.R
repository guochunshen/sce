#'
#'Sends an Email with Status Information
#'
#'@param recipient A valid email address
#'@param subject Subject of the email message
#'@param message Body of the email message
#'
#'
#'@details
#'This is a sample status notification way by sending email from \code{mail} package. 
#'
#'@note
#'To avoid spam, you have a maximum of 20 email notifications per day, that you can send from the same IP-address.
#'
#'@author
#'Dr. Lin Himmelmann <Rmail@linhi.com>
#'


sendmail <- function(recipient, subject="Notification from R",
                    message="Calculation finished!")
{
  password="rmail"
  # Leerteichen kodiert als <
  # Zeilenumbruch "\n" oder "\r" kodiert als >
  pattern = "[^0-9,a-z,A-Z,.,:,;_!,?,%,<,>,=,@,*]"
  recipient = gsub(" "    ,"" , recipient)
  recipient = gsub(pattern,"_", recipient)
  password  = gsub(pattern,"_", password)
  subject   = gsub("[<]"  , "_", subject)
  subject   = gsub(" "    , "<", subject)
  subject   = gsub(pattern, "_", subject)
  message   = gsub("[<>]", "_", message)
  message   = gsub(" "   , "<", message)
  message   = gsub("\n"  , ">", message)
  message   = gsub("\r"  , ">", message)
  message   = gsub(pattern,"_", message)
  path = paste("http://rmail.linhi.de/rmail.php?Email=",recipient,
               "&Passwort=",password,
               "&Betreff=",subject,
               "&Nachricht=",message,
               sep="")
  x = readLines(path)
  return(x)
}
