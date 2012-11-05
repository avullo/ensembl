-- MySQL dump 10.13  Distrib 5.1.61, for redhat-linux-gnu (x86_64)
--
-- Host: mysql-eg-devel-1.ebi.ac.uk    Database: homo_sapiens_core_test_db
-- ------------------------------------------------------
-- Server version	5.1.49-log

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `mapping_session`
--

DROP TABLE IF EXISTS `mapping_session`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `mapping_session` (
  `mapping_session_id` int(11) NOT NULL AUTO_INCREMENT,
  `old_db_name` varchar(80) COLLATE latin1_bin NOT NULL DEFAULT '',
  `new_db_name` varchar(80) COLLATE latin1_bin NOT NULL DEFAULT '',
  `old_release` varchar(5) COLLATE latin1_bin NOT NULL DEFAULT '',
  `new_release` varchar(5) COLLATE latin1_bin NOT NULL DEFAULT '',
  `old_assembly` varchar(20) COLLATE latin1_bin NOT NULL DEFAULT '',
  `new_assembly` varchar(20) COLLATE latin1_bin NOT NULL DEFAULT '',
  `created` datetime NOT NULL,
  PRIMARY KEY (`mapping_session_id`)
) ENGINE=MyISAM AUTO_INCREMENT=4 DEFAULT CHARSET=latin1 COLLATE=latin1_bin;
/*!40101 SET character_set_client = @saved_cs_client */;

/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2012-11-05 10:52:11
