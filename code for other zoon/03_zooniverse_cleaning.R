library(dplyr)
library(tidyr)
basic_text<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/zooniverse_data/limno/text_reducer_basic_LIMN_texts.csv") %>%
  select(subject_id, task, 'data.aligned_text', 'data.number_views', 'data.consensus_score','data.consensus_text')
basic_dropdown<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/zooniverse_data/limno/dropdown_reducer_basic_LIMN_dropdowns.csv") %>%
        select(subject_id, task, 'data.value')
limno_text<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/zooniverse_data/limno/text_reducer_limn_texts_texts.csv") %>%
  select(subject_id, task, 'data.aligned_text', 'data.number_views', 'data.consensus_score','data.consensus_text')
limno_quest<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/zooniverse_data/limno/question_reducer_limn_texts_questions.csv") %>%
  select(subject_id, task, 'data.no.back.side.to.card', 'data.no', 'data.yes')
#add  talk files (e.g. notes that people wrote)
#bluegill 
limno_text_values<- limno_text %>% drop_na('data.number_views')

data_limno<-left_join(basic_dropdown, basic_text, by="subject_id") %>%
          left_join(limno_text_values, by="subject_id") %>%
          left_join(limno_quest, by="subject_id")
#reduced by 800K lines 
write.csv(data_limno, "/Users/katelynking/Desktop/UofM/CHANGES/zooniverse_data/limno/limno_dat.csv", row.names = FALSE)
