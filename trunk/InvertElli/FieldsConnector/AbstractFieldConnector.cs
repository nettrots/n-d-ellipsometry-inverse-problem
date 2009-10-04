using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FieldsConnector
{
    public abstract class AbstractFieldConnector
    {
        private string _name;

        public string Name
        {
            get { return _name; }
            set { _name = value; }
        }

        protected abstract void syncTo();
        protected abstract void syncFrom();


        
    }
 
   

}
/*
 connect(Form.Textbox.Text, Instance.Property,name(string) )
 FormToSearchCnn[name].syncTo()
 FormToSearchCnn[name].syncFrom()
 */