using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FieldsConnector
{
    public abstract class AbstractFieldConnector
    {
        private string _name;
        protected object _field;
        protected object _interface;
      
        public string Name
        {
            get { return _name; }
            set { _name = value; }
        }

        public abstract void syncTo();
        public abstract void syncFrom();

        public bool checkForInterface(object o)
        {
            return o == _interface;
        }

        
    }
 
   

}
/*
 connect(Form.Textbox.Text, Instance.Property,name(string) )
 FormToSearchCnn[name].syncTo()
 FormToSearchCnn[name].syncFrom()
 */